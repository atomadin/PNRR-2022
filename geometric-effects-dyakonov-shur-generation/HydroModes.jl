#=======================================================================
# HydroModes
=======================================================================#
module HydroModes

using LinearAlgebra, OffsetArrays, SparseArrays, Interpolations
using Arpack, ArnoldiMethod, LinearMaps

export Sample, solvesteady, instability

struct ShiftAndInvert{TA,TB,TT}
	A_lu::TA
	B::TB
	temp::TT
end
function (M::ShiftAndInvert)(y,x)
	mul!(M.temp, M.B, x)
	ldiv!(y, M.A_lu, M.temp)
end
function construct_linear_map(A,B)
	a = ShiftAndInvert(lu(A),B,Vector{eltype(A)}(undef, size(A,1)))
	LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
end

use_arpack=true

#=======================================================================
# Constants
=======================================================================#
const δx, δy = CartesianIndex(1,0), CartesianIndex(0,1)

#=======================================================================
# Types
=======================================================================#

#=======================================================================
A contact is represented as a `Tuple` of:
- `Int`: index of the boundary;
- `Int`, `Int`: index of the first and last cell;
- `Char`: type, potential 'V', current 'I';
- `T`: value of the fixed quantity.
=======================================================================#
Contact{T<:Number} = Tuple{Int, Int, Int, Char, T}

#=======================================================================
Boundary direction: sample to the left.
A boundary is represented as a `Tuple` of:
- `Int`, `Int`: first and last $x$ coordinate;
- `Int`, `Int`: first and last $y$ coordinate;
- `Int`, `Int`: index of the first and last contact.
=======================================================================#
Boundary = NTuple{6, Int}

#=======================================================================
A sample is represented as struct consisting of:
- `Nx::Int`, `Ny::Int`: divisions along $x$, $y$;
- `L::T`, `W::T`: sample length ($x$) and width ($y$);
- `s::T`: plasma-wave velocity;
- `τ::T`: momentum relaxation time;
- `ν::T`: viscosity;
- `lb::T`: slip length;
- `boundaries::Vector{Boundary}`: list of boundaries;
- `contacts::Vector{Contact{T}}`: list of boundaries;
=======================================================================#
Base.Base.@kwdef struct Sample{T<:Number}
	Nx::Int; Ny::Int
	L::T = 1; W::T
	s::T = 1
	τ::T
	ν::T
	lb::T = Inf
	boundaries::Vector{Boundary}
	contacts::Vector{Contact{T}}
end

#=======================================================================
# Methods
=======================================================================#

#=======================================================================
## Variables and indices
=======================================================================#
function indices(sample::Sample)
	Nx, Ny = sample.Nx, sample.Ny
	return OffsetArray(LinearIndices((Nx+2, Ny+2, 3)), -1, -1, 0)
end
function variables(sample::Sample)
	return eachslice(indices(sample), dims=3)
end

function dynamical(sample::Sample)
	Nx, Ny = sample.Nx, sample.Ny
	return (
		CartesianIndices((1:Nx, 1:Ny)),
		CartesianIndices((1:Nx-1, 1:Ny)),
		CartesianIndices((1:Nx, 1:Ny-1))
	)
end
function indicesdyn(sample::Sample)
	n, u, v = variables(sample)
	rn, ru, rv = dynamical(sample)
	return [ vec(n[rn]); vec(u[ru]); vec(v[rv]) ]
end

function pusheq!(eqs, vjs::Pair...)
	for vj in vjs
		for (v, j) in zip(vj[1], vj[2])
			push!(eqs.I, eqs.eq[1])
			push!(eqs.J, j)
			push!(eqs.V, v)
		end
	end
end

function dyneqs(
	sample::Sample{T}, x0::AbstractArray{T,3}
	) where {T}
	# parameters
	Nx, Ny, L, W, s, τ, ν = getproperty.((sample,),
		(:Nx, :Ny, :L, :W, :s, :τ, :ν))
	Dx, Dy = Nx / L, Ny / W

	# variables
	n, u, v = eachslice(x0, dims=3)
	rn, ru, rv = dynamical(sample)

	# equations
	eqs = T[]
	sizehint!(eqs, length(rn) + length(ru) + length(rv))
	
	for r in rn
		push!(eqs,
		Dx/2 * ((n[r] + n[r+δx]) * u[r] - (n[r-δx] + n[r]) * u[r-δx]) +
		Dy/2 * ((n[r] + n[r+δy]) * v[r] - (n[r-δy] + n[r]) * v[r-δy])
		)
	end

	for r in ru
		push!(eqs,
		u[r] * Dx/2 * (u[r+δx] - u[r-δx]) +
		(v[r] + v[r+δx] + v[r+δx-δy] + v[r-δy])/4 *
			Dy/2 * (u[r+δy] - u[r-δy]) +
		s^2 * Dx * (n[r+δx] - n[r]) +
		1/τ * u[r] -
		ν * Dx^2 * (u[r-δx] - 2u[r] + u[r+δx]) -
		ν * Dy^2 * (u[r-δy] - 2u[r] + u[r+δy])
		)
	end

	for r in rv
		push!(eqs,
		(u[r] + u[r-δx] + u[r-δx+δy] + u[r+δy])/4 *
			Dx/2 * (v[r+δx] - v[r-δx]) +
		v[r] * Dy/2 * (v[r+δy] - v[r-δy]) +
		s^2 * Dy * (n[r+δy] - n[r]) +
		1/τ * v[r] -
		ν * Dx^2 * (v[r-δx] - 2v[r] + v[r+δx]) -
		ν * Dy^2 * (v[r-δy] - 2v[r] + v[r+δy])
		)
	end
	
	return eqs
end

function dynjac(
	sample::Sample{T}, x0::AbstractArray{T,3}
	) where {T}
	# parameters
	Nx, Ny, L, W, s, τ, ν = getproperty.((sample,),
		(:Nx, :Ny, :L, :W, :s, :τ, :ν))
	Dx, Dy = Nx / L, Ny / W

	# variables
	n0, u0, v0 = eachslice(x0, dims=3)
	n, u, v = variables(sample)
	rn, ru, rv = dynamical(sample)

	# equations
	eqs = (eq=[0], I = Int[], J = Int[], V = T[])
	sizehint!.((eqs.I, eqs.J, eqs.V), (10+18+18) * Nx * Ny)
	
	for r in rn
		eqs.eq[1] += 1
		pusheq!(eqs,
		Dx/2 .* (n0[r]+n0[r+δx], -n0[r-δx]-n0[r]) => (u[r], u[r-δx]),
		Dx/2 .* (u0[r]-u0[r-δx], u0[r], -u0[r-δx]) => (n[r], n[r+δx], n[r-δx]),
		Dy/2 .* (n0[r]+n0[r+δy], -n0[r-δy]-n0[r]) => (v[r], v[r-δy]),
		Dy/2 .* (v0[r]-v0[r-δy], v0[r], -v0[r-δy]) => (n[r], n[r+δy], n[r-δy]),
		)
	end

	for r in ru
		eqs.eq[1] += 1
		pusheq!(eqs,
		u0[r] * Dx/2 .* (1, -1) => (u[r+δx], u[r-δx]),
		Dx/2 * (u0[r+δx] - u0[r-δx]) => u[r],
		(v0[r] + v0[r+δx] + v0[r+δx-δy] + v0[r-δy])/4 * Dy/2 .* (1, -1) =>
			(u[r+δy], u[r-δy]),
		Dy/2 * (u0[r+δy] - u0[r-δy]) / 4 .* (1, 1, 1, 1) =>
			(v[r], v[r+δx], v[r+δx-δy], v[r-δy]),
		s^2 * Dx .* (1, -1) => (n[r+δx], n[r]),
		1/τ => u[r],
		-ν * Dx^2 .* (1, -2, 1) => (u[r-δx], u[r], u[r+δx]),
		-ν * Dy^2 .* (1, -2, 1) => (u[r-δy], u[r], u[r+δy]),
		)
	end

	for r in rv
		eqs.eq[1] += 1
		pusheq!(eqs,
		(u0[r] + u0[r-δx] + u0[r-δx+δy] + u0[r+δy])/4 * Dx/2 .* (1, -1) =>
			(v[r+δx], v[r-δx]),
		Dx/2 * (v0[r+δx] - v0[r-δx]) / 4 .* (1, 1, 1, 1) =>
			(u[r], u[r-δx], u[r-δx+δy], u[r+δy]),
		v0[r] * Dy/2 .* (1, -1) => (v[r+δy], v[r-δy]),
		Dy/2 * (v0[r+δy] - v0[r-δy]) => v[r],
		s^2 * Dy .* (1, -1) => (n[r+δy], n[r]),
		1/τ => v[r],
		-ν * Dx^2 .* (1, -2, 1) => (v[r-δx], v[r], v[r+δx]),
		-ν * Dy^2 .* (1, -2, 1) => (v[r-δy], v[r], v[r+δy]),
		)
	end
	
	return sparse(eqs.I, eqs.J, eqs.V, eqs.eq[1], length(x0))
end

function bcseqs(
	sample::Sample{T}, x0::AbstractArray{T,3}
	) where {T}
	# parameters
	Nx, Ny, L, W, lb = getproperty.((sample,),
		(:Nx, :Ny, :L, :W, :lb))
	Dx, Dy = Nx / L, Ny / W

	# variables
	n, u, v = eachslice(x0, dims=3)

	# equations
	eqs = T[]
	sizehint!(eqs, 3 * (Nx + Ny))

	l1, l2 = lb==Inf ? (1, 0) : (lb, 1) 

	for (x1, x2, y1, y2, c1, c2) in sample.boundaries
		# boundary coordinates
		r1, r2 = extrema((CartesianIndex(x1,y1), CartesianIndex(x2,y2)))
		rb = r1:r2
		# tangential and normal directions without sign
		δn, δt = x1==x2 ? (δx, δy) : (δy, δx)
		# differential with sign
		Dt = sign(x2-x1) * Dx + sign(y2-y1) * Dy
		Dn = sign(y2-y1) * Dx - sign(x2-x1) * Dy
		# tangential and normal velocities and indices without sign
		vn, vt = x1==x2 ? (u, v) : (v, u)
		
		# surface boundary condition
		for r in @view rb[begin:end-1]
			push!(eqs,
			Dx * (v[r+δx] - v[r]) +
			Dy * (u[r+δy] - u[r])
			#l1 * sign(Dn) * Dt * (vn[r+δt] - vn[r]) +
			#l1 * sign(Dt) * Dn * (vt[r+δn] - vt[r]) +
			#l2 * sign(Dt) /  2 * (vt[r+δn] + vt[r])
			)
		end

		# edge mask
		edge = trues(size(rb))
		# contacts
		for (_, i1, i2, type, value) in @view sample.contacts[c1:c2]
			rc = @view rb[i1:i2]
			edge[i1:i2] .= false
			if type == 'V'
				for r in rc
					push!(eqs, (n[r] + n[r+δn])/2 - value)
				end
			elseif type == 'I'
				for (ri, rj) in zip(rc[begin:end-1], rc[begin+1:end])
					push!(eqs,
					(n[ri] + n[ri+δn])/2 -
					(n[rj] + n[rj+δn])/2
					)
				end
				j0 = zero(T)
				for r in rc
					j0 += (n[r] + n[r+δn])/2 * vn[r]
				end
				push!(eqs, j0 - value * length(rc))
			end
		end
		for r in rb[edge]
			push!(eqs, vn[r])
		end

	end
	
	rn =  dynamical(sample)[1]
	for r in rn
		if !(r-δx in rn) || !(r+δx in rn)
			push!(eqs, 2n[r] - n[r-δx] - n[r+δx])
		end
		if !(r-δy in rn) || !(r+δy in rn)
			push!(eqs, 2n[r] - n[r-δy] - n[r+δy])
		end
	end
	
	push!(eqs,
	n[0,0], n[Nx+1,0], n[Nx+1,Ny+1], n[0,Ny+1],
	u[0,0], u[Nx,0], u[0,Ny+1], u[Nx,Ny+1],
	v[0,0], v[0,Ny], v[Nx+1,0], v[Nx+1,Ny])
	for i in u[Nx+1,:]
		push!(eqs, i)
	end
	for i in v[:,Ny+1]
		push!(eqs, i)
	end
	
	return eqs
	
end

function bcsjac(
	sample::Sample{T}, x0::AbstractArray{T,3}
	) where {T}
	# parameters
	Nx, Ny, L, W, lb = getproperty.((sample,),
		(:Nx, :Ny, :L, :W, :lb))
	Dx, Dy = Nx / L, Ny / W

	# variables
	n0, u0, v0 = eachslice(x0, dims=3)
	n, u, v = variables(sample)

	# equations
	eqs = (eq=[0], I = Int[], J = Int[], V = T[])
	sizehint!.((eqs.I, eqs.J, eqs.V), 3 * (Nx + Ny))
	
	l1, l2 = lb==Inf ? (1, 0) : (lb, 1) 

	for (x1, x2, y1, y2, c1, c2) in sample.boundaries
		# boundary coordinates
		r1, r2 = extrema((CartesianIndex(x1,y1), CartesianIndex(x2,y2)))
		rb = r1:r2
		# tangential and normal directions without sign
		δn, δt = x1==x2 ? (δx, δy) : (δy, δx)
		# differential with sign
		Dt = sign(x2-x1) * Dx + sign(y2-y1) * Dy
		Dn = sign(y2-y1) * Dx - sign(x2-x1) * Dy
		# tangential and normal velocities and indices without sign
		vn, vt = x1==x2 ? (u, v) : (v, u)
		v0n, v0t = x1==x2 ? (u0, v0) : (v0, u0)
		
		# surface boundary condition
		for r in @view rb[begin:end-1]
			eqs.eq[1] += 1
			pusheq!(eqs,
			Dx .* (1, -1) => (v[r+δx], v[r]),
			Dy .* (1, -1) => (u[r+δy], u[r]),
			#l1 * sign(Dn) * Dt .* (1, -1) => (vn[r+δt], vn[r]),
			#l1 * sign(Dt) * Dn .* (1, -1) => (vt[r+δn], vt[r]),
			#l2 * sign(Dt) /  2 .* (1,  1) => (vt[r+δn], vt[r]),
			)
		end

		# edge mask
		edge = trues(size(rb))
		# contacts
		for (_, i1, i2, type, value) in @view sample.contacts[c1:c2]
			rc = @view rb[i1:i2]
			edge[i1:i2] .= false
			if type == 'V'
				for r in rc
					eqs.eq[1] += 1
					pusheq!(eqs, (1/2, 1/2) => (n[r], n[r+δn]))
				end
			elseif type == 'I'
				for (ri, rj) in zip(rc[begin:end-1], rc[begin+1:end])
					eqs.eq[1] += 1
					pusheq!(eqs,
					( 1/2,  1/2) => (n[ri], n[ri+δn]),
					(-1/2, -1/2) => (n[rj], n[rj+δn]),
					)
				end
				eqs.eq[1] += 1
				for r in rc
					pusheq!(eqs,
					v0n[r] .* (1/2, 1/2) => (n[r], n[r+δn]),
					(n0[r] + n0[r+δn])/2 => vn[r],
					)
				end
			end
		end
		for r in rb[edge]
			eqs.eq[1] += 1
			pusheq!(eqs, 1 => vn[r])
		end

	end
	
	rn = dynamical(sample)[1]
	for r in rn
		if !(r-δx in rn) || !(r+δx in rn)
			eqs.eq[1] += 1
			pusheq!(eqs, (2, -1, -1) => (n[r], n[r-δx], n[r+δx]))
		end
		if !(r-δy in rn) || !(r+δy in rn)
			eqs.eq[1] += 1
			pusheq!(eqs, (2, -1, -1) => (n[r], n[r-δy], n[r+δy]))
		end
	end
	
	for j in (
		n[0,0], n[Nx+1,0], n[Nx+1,Ny+1], n[0,Ny+1],
		u[0,0], u[Nx,0], u[0,Ny+1], u[Nx,Ny+1],
		v[0,0], v[0,Ny], v[Nx+1,0], v[Nx+1,Ny],
		)
		eqs.eq[1] += 1
		pusheq!(eqs, 1 => j)
	end
	for j in u[Nx+1,:]
		eqs.eq[1] += 1
		pusheq!(eqs, 1 => j)
	end
	for j in v[:,Ny+1]
		eqs.eq[1] += 1
		pusheq!(eqs, 1 => j)
	end
	
	return sparse(eqs.I, eqs.J, eqs.V, eqs.eq[1], length(x0))
end

tolerance(x1, x0) = sqrt(sum(abs2, x1) / sum(abs2, x0))

function solvesteady(sample::Sample{T}; maxiter=0, tol=1e-3) where {T}
	Nx, Ny = sample.Nx, sample.Ny
	γ = 1e0

	# initial state
	x0 = similar(indices(sample), T)
	n0s, u0s, v0s = one(T), zero(T), zero(T)
	x0[:,:,1] .= n0s; x0[:,:,2] .= u0s; x0[:,:,3] .= v0s
	#x0 = reshape(randn(T, size(indices(sample))), axes(indices(sample)))

	# 0-th iteration
	E0 = vcat(dyneqs(sample, x0), bcseqs(sample, x0))
	J0 = vcat(dynjac(sample, x0), bcsjac(sample, x0))
	# LU decomposition
	F = lu(J0)
	# qr decomposition maybe faster for rank(qr(AB))
	if rank(F.U) != size(J0, 1)
		@error "Inconsistent boundary conditions (non-full rank)"
	end
	# condition for a fixed point: |det(I + A⁻¹)| = |det(I + A) / det(A)| < 1
	# approximation: det(I + ϵA) = 1 + ϵ tr(A) + O(ϵ²)
	if abs((1 + tr(J0)) / det(F)) >= 1
		@warn "The iteration might not have a fixed point"
	end
	# update
	vec(x0) .-= γ .* (F \ E0)
	# tolerance
	E0 = vcat(dyneqs(sample, x0), bcseqs(sample, x0))
	t = sqrt(sum(abs2, E0)) / length(E0)
	#t = sqrt(sum(abs2, E0)) / (3Nx+3Ny)
	@info "iteration 0" tolerance=t

	# iterations
	for i in 1:maxiter
		J0 = vcat(dynjac(sample, x0), bcsjac(sample, x0))
		# use previous factorization
		F = lu!(F, J0)
		# fixed point
		if abs((1 + tr(J0)) / det(F)) >= 1
			@warn "The iteration might not have a fixed point"
		end
		# update
		γ /= 2
		vec(x0) .-= γ .* (F \ E0)
		# tolerance and warning
		E0 = vcat(dyneqs(sample, x0), bcseqs(sample, x0))
		t = sqrt(sum(abs2, E0)) / length(E0)
		#t = sqrt(sum(abs2, E0)) / (3Nx+3Ny)
		@info "iteration $i" tolerance=t
		if t <= tol break
		elseif i == maxiter
			@warn "Maximum iteration number reached" maxiter=maxiter tolerance=t
		end
	end

	return x0
end

function q(F)
	n, m = size(F.factors)
	q = sparse(one(eltype(F.τ))*I, n, n)
	for i in 1:m
		q -= F.τ[i] .* (q * F.factors[:,i] * F.factors[:,i]')
	end
	return permute!(q, F.rpivinv, 1:n)
end

function modes(sample::Sample{T}, steady::NTuple{3, AbstractArray{T,3}};
		nev=10, tol=1e-5, maxiter=300, sigma=nothing) where {T}
	# linearization
	A, _ = lineareqs(sample, steady)
	B, _ = linearbcs(sample, steady)
	# number of variables
	Nvar = size(A, 2)
	# number of boundary conditions
	Nbcs = size(B, 1)
	# indices of dynamical variables
	D = indicesdyn(sample)
	# number of dynamical variables
	Ndyn = length(D)

	@time "QR" begin
		qrB = qr(sparse(B'))
		Q = q(qrB)
	end

	A0 = (A * Q)[:, Nbcs+1:end]
	Q0 = Q[D, Nbcs+1:end]
	# immersion
	P = sparse(D, 1:Ndyn, ones(T, Ndyn), Nvar, Ndyn)

	if use_arpack
	@time "F" F = eigs(complex(Q0), complex(A0),
		nev=nev, which=:LM, tol=tol, maxiter=maxiter, sigma=sigma, ritzvec=true)
	else
	@time "F" begin
		# Target the largest eigenvalues of the inverted problem
		decomp,  = partialschur(construct_linear_map(complex(A0), complex(Q0)),
			nev=nev, tol=tol, restarts=maxiter, which=LM())
		F = partialeigen(decomp)
	end
	end

	return -im ./ F[1], P * Q0 * F[2]
end

function instability(sample, steady; kwargs...)
	ω0, M0 = modes(sample, steady; kwargs...)
	unstable = @. imag(ω0) >= 0 && real(ω0) >= 0
	ω = ω0[unstable]
	M = [ reshape(c, axes(indices(sample)))
		for c in eachcol(@view M0[:, unstable]) ]
	return ω, M
end

function ipinsta(sample, steady; kwargs...)
	ω, x = instability(sample, steady; kwargs...)
	return ω, interpolation.((sample,), x)
end

end
