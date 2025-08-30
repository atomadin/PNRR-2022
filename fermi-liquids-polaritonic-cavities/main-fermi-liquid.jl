using OffsetArrays

abstract type Dielectric{T} end
struct UniaxialDielectric{T,N} <: Dielectric{T}
	εxInf::T; sx::NTuple{N,T}; ωxT::NTuple{N,T}; γx::NTuple{N,T};
	εzInf::T; sz::NTuple{N,T}; ωzT::NTuple{N,T}; γz::NTuple{N,T};
end

struct Device{T,N}
	dielectric::UniaxialDielectric{T,N}
	h1::T; h2::T
end

# units:
#  frequency: meV
#  length: nm
#  density: 1e12 cm^-2

const hBN = UniaxialDielectric{Float64,1}(
	5.1, 0.124 .* (1994,), 0.124 .* (1394,), 0.124 .* (1.8,),
	2.5, 0.124 .* (495,),  0.124 .* (785,),  0.124 .* (1,),
)

const Bi2Se3 = UniaxialDielectric{Float64,2}(
	29.0, 0.124 .* (704, 55),  0.124 .* (64,  125), 0.124 .* (3.5, 3.5), 
	17.4, 0.124 .* (283, 156), 0.124 .* (135, 154), 0.124 .* (3.5, 3.5),
)

const Bi2Te3 = UniaxialDielectric{Float64,2}(
	85.0, 0.124 .* (667, 181), 0.124 .* (50, 95),  0.124 .* (3.5, 10),
	50.0, 0.124 .* (314, 353), 0.124 .* (94, 120), 0.124 .* (3.5, 15),
)

const Nf = 4e0
const αee = 300e0 / 137e0
const ħvF = 6.582119569e-13 * 1e15

kFωF(n0) = sqrt(1e-2 * π * n0) .* (1e0, ħvF)

function ε(ω::T, m::UniaxialDielectric{T,N}) where {T,N}
	εx = m.εxInf; εz = m.εzInf
	for i in 1:N
		εx += m.sx[i]^2 / (m.ωxT[i]^2 + ω^2 + m.γx[i] * ω)
		εz += m.sz[i]^2 / (m.ωzT[i]^2 + ω^2 + m.γz[i] * ω)
	end
	return εx, εz
end

function χ0(q::T, ω::T) where {T}
	C1 = q^2 / sqrt(ω^2 + q^2)
	C2 = (2 + im * ω) / q
	return -Nf / (2π) * (1 + C1 * (π/8 - 1/4 * real(asin(C2) + C2 * sqrt(1 - C2^2))))
end

# NOTE: the F... functions take as arguments q and ω in physical units

Fvacuum(q::T, ω::T, d::Device{T}) where {T} = one(T)

function Fdynamic(q::T, ω::T, d::Device{T}) where {T}
	εx, εz = ε(ω, d.dielectric)
	sqrtεxεz = sqrt(εx / εz)
	C0 = 2 / εz / sqrtεxεz
	C1 = exp(-2 * sqrtεxεz * q * d.h1)
	C2 = exp(-2 * sqrtεxεz * q * d.h2)
	return C0 / ((1 + C1) / (1 - C1) + (1 + C2) / (1 - C2))
end

Fstatic(q, ω, d) = Fdynamic(q, zero(ω), d)

function tanhsinh(N, h)
	π_2 = π / 2
	x, w = zeros(-N:N), zeros(-N:N)
	for i in -N:N
		x[i] = tanh(π_2 * sinh(i * h))
		w[i] = h * π_2 * cosh(i * h) * (1 - x[i]^2)
	end
	return x, w
end

function expsinh(N, h)
	π_2 = π / 2
	x, w = zeros(-N:N), zeros(-N:N)
	for i in -N:N
		x[i] = exp(π_2 * sinh(i * h))
		w[i] = h * π_2 * cosh(i * h) * x[i]
	end
	return x, w
end

# compute derivatives of the self-energy
function ∂ωΣ∂kΣ(s1, n0::T, F::Function, d::Device{T}) where {T}
	kF, ωF = kFωF(n0)
	Λ = sqrt(4π/3/sqrt(3)) / (0.142 * kF)
	Nω, Nq, Nφ = 50, 50, 50
	xω, wω = expsinh(Nω, 3e0 / Nω)
	xq, wq = tanhsinh(Nq, 3e0 / Nq)
	xφ, wφ = tanhsinh(Nφ, 2e0 / Nφ)
	@inline function ∫_res(q)
		Fq0 = F(q * kF, zero(ωF), d)
		εq0 = 1 + Nf * αee * Fq0 / q
		return Fq0 / q / εq0 * (1 + s1 * (1 - q^2 / 2)) / sqrt(1 - q^2 / 4)
	end
	function ∫dq_res()
		I1 = zero(T)
		for i in eachindex(xq)
			I1 += wq[i] * ∫_res(1 + xq[i])
		end
		return I1
	end
	@inline function ∫_line(q, φ, ω)
		Fqω = F(q * kF, ω * ωF, d)
		εqω = 1 - 2π * αee / q * Fqω * χ0(q, ω)
		kq = sqrt(1 + q^2 - 2 * q * cos(φ))
		cosθ = (1 - q * cos(φ)) / kq
		∂kcosθ = q^2 * (1 - cos(φ)^2) / kq^3
		ξ(s2) = s2 * kq - 1
		I1(s2) = Fqω / εqω * (1 + s1 * s2 * cosθ) * (ω^2 - ξ(s2)^2) / (ω^2 + ξ(s2)^2)^2
		J1(s2) = Fqω / εqω * (s1 * s2 * ∂kcosθ * ξ(s2) / (ω^2 + ξ(s2)^2) +
			s2 * (1 + s1 * s2 * cosθ) * cosθ * (ω^2 - ξ(s2)^2) / (ω^2 + ξ(s2)^2)^2)
		return (I1(1) + I1(-1), J1(1) + J1(-1))
	end
	@inline function ∫dω_line(q, φ)
		I1 = J1 = zero(T)
		for i in eachindex(xω)
			i1, j1 = ∫_line(q, φ, xω[i])
			I1 += wω[i] * i1
			J1 += wω[i] * j1
		end
		return I1, J1
	end
	@inline function ∫dqdω_line(φ)
		cosφ = cos(φ)
		if 0 < cosφ < Λ/2
			C1 = Λ/2 - cosφ
			C2 = Λ/2 + cosφ
			I1 = J1 = I2 = J2 = zero(T)
			for i in eachindex(xq)
				i1, j1 = ∫dω_line(cosφ * (1 + xq[i]), φ)
				I1 += cosφ * wq[i] * i1
				J1 += cosφ * wq[i] * j1
				i2, j2 = ∫dω_line(C2 + C1 * xq[i], φ)
				I2 += C1 * wq[i] * i2
				J2 += C1 * wq[i] * j2
			end
			return I1 + I2, J1 + J2
		else
			I1 = J1 = zero(T)
			for i in eachindex(xq)
				i1, j1 = ∫dω_line(Λ/2 * (1 + xq[i]), φ)
				I1 += Λ/2 * wq[i] * i1
				J1 += Λ/2 * wq[i] * j1
			end
			return I1, J1
		end
	end
	function ∫dφdqdω_line()
		π_4 = π / 4
		I1 = I2 = J1 = J2 = zero(T)
		for i in eachindex(xφ)
			i1, j1 = ∫dqdω_line(π_4 * (1 + xφ[i]))
			I1 += π_4 * wφ[i] * i1
			J1 += π_4 * wφ[i] * j1
			i2, j2 = ∫dqdω_line(π_4 * (3 + xφ[i]))
			I2 += π_4 * wφ[i] * i2
			J2 += π_4 * wφ[i] * j2
		end
		return 2 * (I1 + I2), 2 * (J1 + J2)
	end
	∂ω_line, ∂k_line = ∫dφdqdω_line()
	∂ω_res = ∫dq_res()
	return (αee / 2π) .* (∂ω_res - ∂ω_line / 2π, ∂k_line / 2π)
end

# renormalized Fermi velocity
function vFDyson(∂ωΣ, ∂kΣ)
	return (1 + ∂kΣ) / (1 - ∂ωΣ)
end

