# constants
const me = 9.1e-28
const e = 4.8e-10
const ħ = 1e-27
# GaAs
const mb = 0.063me
# graphene
const vF = 1e8

struct DSDevice
	L # channel length (cm)
	Cb; Ct; Cq; C1; C2 # capacitance bottom, top, quantum (cm^-1)
	s; s1; s2 # plasma velocities (cm/s)
	ω1; ω2 # plasma frequencies (Hz)
	n0 # density (cm^-2)
	τ # momentum relaxation time (s)
	ν # viscosity (cm^2/s)
	η # parab: 0, conic: 1/2
	function DSDevice(L, Cb, Ct, Cq, m, n0, τ, ν, η)
		Vb = e * n0 / Cb
		s = sqrt(e * Vb / m)
		C1 = (Cb^-1 + Cq^-1)^-1
		C2 = ((Cb+Ct)^-1 + Cq^-1)^-1
		s1 = s * sqrt(Cb / C1); ω1 = π * s1 / 2L
		s2 = s * sqrt(Cb / C2); ω2 = π * s2 / 2L
		new(L, Cb, Ct, Cq, C1, C2, s, s1, s2, ω1, ω2, n0, τ, ν, η)
	end
end

"""
    semicond_device(L, Cb, Ct, n0, τ, ν)

Gallium arsenide Dyakonov-Shur photodetection device as a `DSDevice`.

# Arguments
- `L`: channel length (cm).
- `Cb`: bottom capacitance `εb/(4π * hb)` (cm⁻¹).
- `Ct`: top capacitance `εt/(4π * ht)` (cm⁻¹).
- `n0`: 2D electron gas average density (cm⁻²).
- `τ`: momentum relaxation time (s).
- `ν`: 2D electron gas viscosity (cm² s⁻¹).
"""
function semicond_device(L, Cb, Ct, n0, τ, ν)
	Cq = mb * e^2 / ħ^2 / π
	DSDevice(L, Cb, Ct, Cq, mb, n0, τ, ν, 0)
end

"""
    graphene_device(L, Cb, Ct, n0, τ, ν)

Single-layer graphene Dyakonov-Shur photodetection device as a `DSDevice`.

# Arguments
- `L`: channel length (cm).
- `Cb`: bottom capacitance `εb/(4π * hb)` (cm⁻¹).
- `Ct`: top capacitance `εt/(4π * ht)` (cm⁻¹).
- `n0`: 2D electron gas average density (cm⁻²).
- `τ`: momentum relaxation time (s).
- `ν`: 2D electron gas viscosity (cm² s⁻¹).
"""
function graphene_device(L, Cb, Ct, n0, τ, ν)
	mc = ħ * sqrt(π * n0) / vF
	Cq = 2 * e^2 * n0 / ħ / vF / sqrt(π * n0)
	DSDevice(L, Cb, Ct, Cq, mc, n0, τ, ν, 1/2)
end

K1(D::DSDevice, ω) = sqrt(ω * (ω + im/D.τ) / (D.s1^2 - im * D.ν * ω))
K2(D::DSDevice, ω) = sqrt(ω * (ω + im/D.τ) / (D.s2^2 - im * D.ν * ω))

"""
    coeff(D::DSDevice, ω, l::Vector)

Compute the coefficients of the first-order solution for the device `D`, with regions of length `l`, at the frequency `ω`.
"""
function coeff(D::DSDevice, ω, l::Vector)
	n = length(l) # 2N, for N gates
	Σl = sum(l)
	@assert iseven(n)
	@assert Σl ≤ 1
	@assert all(l .≥ 0)
	kb = K1(D, ω)
	kt = K2(D, ω)
	ktkb = kt / kb
	ua = D.Ct / (D.Cb + D.Ct)
	c1 = D.Cb / D.C1 - im * ω * D.ν / D.s^2
	c2 = D.Cb / D.C2 - im * ω * D.ν / D.s^2
	M = zeros(typeof(kb), 2n+2, 2n+2)
	m = zeros(typeof(kb), 2n+2)
	for j in 1:2:n
		ikb = im * kb * D.L * l[j]
		ikt = im * kt * D.L * l[j+1]
		# v_1(j), v_1(j+1)
		@views M[2j-1, 2j-1:2j+2] .= ktkb * exp(ikb), -ktkb * exp(-ikb), -1, 1
		# v_1(j+1), v_1(j+2)
		@views M[2j  , 2j+1:2j+4] .= exp(ikt), -exp(-ikt), -ktkb, ktkb
		# n_1(j), n_1(j+1)
		m[2j+1] = -ua
		@views M[2j+1, 2j-1:2j+2] .= c1 * exp(ikb), c1 * exp(-ikb), -c2, -c2
		# n_1(j+1), n_1(j+2)
		m[2j+2] =  ua
		@views M[2j+2, 2j+1:2j+4] .= c2 * exp(ikt), c2 * exp(-ikt), -c1, -c1
	end
	ikb = im * kb * D.L * (1 - Σl)
	@views M[end-1, 1:2] .= 1, 1
	@views M[end, end-1:end] .= exp(ikb), -exp(-ikb)
	return M \ m
end

function ab(j, D::DSDevice, ω)
	s = isodd(j) ? D.s1 : D.s2
	η, ν, τ = D.η, D.ν, D.τ
	sω2ντ = s^2 + ω^2 * ν * τ
	a = 2 - η - ν / τ * (
	(1 + ω^2 * τ^2) / sω2ντ + 2sω2ντ / (s^4 + ω^2 * ν^2))
	b = 1 - η - η/2 * (vF / s)^2 - ω * ν * (ω * τ - ω * ν / s^2) / sω2ντ
	return a, b
end

"""
    ΔU(D::DSDevice, ω, l::Vector)

Compute the photoresponse for the device `D`, with regions of length `l`, at the frequency `ω`.

The top gates divide the channel in regions.
The odd regions are single-gated, the even ones are dual-gated.
The vector `l` represents the length of each region in units of `D.L`.
The last region is not included in `l` and has length `1 - sum(l)`.
"""
function ΔU(D::DSDevice, ω, l::Vector)
	A = coeff(D, ω, l)
	n = length(l)
	Σl = sum(l)
	δφ = zero(ω)
	for j in 1:2:n # odd
		k = K1(D, ω)
		ik = im * k * D.L * l[j]
		a, b = ab(j, D, ω)
		δφ += a * abs2(ω / k / D.s1) * (
		abs2(A[2j-1] - A[2j]) - abs2(A[2j-1] * exp(ik) - A[2j] * exp(-ik)))
		δφ += b * (
		abs2(A[2j-1] + A[2j]) - abs2(A[2j-1] * exp(ik) + A[2j] * exp(-ik)))
	end
	for j in 2:2:n # even
		k = K2(D, ω)
		ik = im * k * D.L * l[j]
		a, b = ab(j, D, ω)
		δφ += a * abs2(ω / k / D.s1) * (
		abs2(A[2j-1] - A[2j]) - abs2(A[2j-1] * exp(ik) - A[2j] * exp(-ik)))
		δφ += b * (D.s2 / D.s1)^2 * (
		abs2(A[2j-1] + A[2j]) - abs2(A[2j-1] * exp(ik) + A[2j] * exp(-ik)))
	end
	# last
	j = n + 1
	k = K1(D, ω)
	ik = im * k * D.L * (1 - Σl)
	a, b = ab(j, D, ω)
	δφ += a * abs2(ω / k / D.s1) * (
	abs2(A[2j-1] - A[2j]) - abs2(A[2j-1] * exp(ik) - A[2j] * exp(-ik)))
	δφ += b * (
	abs2(A[2j-1] + A[2j]) - abs2(A[2j-1] * exp(ik) + A[2j] * exp(-ik)))
	return δφ
end
