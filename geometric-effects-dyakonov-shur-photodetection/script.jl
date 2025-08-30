include("DyakonovShur.jl")
using Printf
# NOTE: cgs units, (tension measured as 1 statV ≈ 300 V)

# geometry
hb = .100e-4 # 100 nm
ht = .100e-4
Cb1 = 1e1 / (4π * hb)
Ct1 = 1e1 / (4π * ht)
Cb2 = 3e0 / (4π * hb)
Ct2 = 3e0 / (4π * ht)
L = 10e-4

K1(D::DSDevice, ω) = sqrt(ω * (ω + im/D.τ) / (D.s1^2 - im * D.ν * ω))

K2(D::DSDevice, ω) = sqrt(ω * (ω + im/D.τ) / (D.s2^2 - im * D.ν * ω))

function ΔUa(D::DSDevice, ω, x0, a)
	k1 = 2D.L * real(K1(D, ω))
	k2 = 2D.L * imag(K1(D, ω))
	cs0 = cos(k1 * x0)
	cs1 = cos(k1 * (1-x0))
	ch0 = cosh(k2 * x0)
	ch1 = cosh(k2 * (1-x0))
	return (D.Ct * D.C2 / D.Cb / (D.Cb + D.Ct) * ω * a * D.L / D.s1)^2 * (
		(2 - D.η) * (ch1 + cs1)
		-(1 - D.η - D.η/2 * (vF / D.s2)^2) *(cs0 * ch1 + cs1 * ch0)
		-(1 - D.η - D.η/2 * (vF / D.s1)^2)
		* sqrt(1 + (ω * D.τ)^-2) * (ch0 - cs0)
	) / (cosh(k2) + cos(k1))
end

file = open("lnu", "w")
d1 = graphene_device(L, Cb2, Ct2, 1e12, 2e-12, 0e0)
d2 = graphene_device(L, Cb2, Ct2, 1e12, 2e-12, 5e2)
for ω in range(1e-3, 6, length=512)
	@printf(file, "%.6g\t%.6g\t%.6g\n", ω,
		ΔU(d1, ω*d1.ω1, [.4, .2]),
		ΔU(d2, ω*d2.ω1, [.4, .2]))
end
close(file)

file = open("papprox", "w")
d = semicond_device(L, Cb1, Ct1,1e11, 2e-11, 1e1)
x0 = 1/2
a = 1/10
for ω in range(1e-3, 6, length=512)
	@printf(file, "%.6g\t%.6g\t%.6g\n", ω,
		ΔU(d, ω*d.ω1, [x0-a/2, a]),
		ΔUa(d, ω*d.ω1, x0, a))
end
close(file)

file = open("ptau", "w")
d0 = semicond_device(L, Cb1, Ct1, 1e11, 5e-11, 1e2)
d = semicond_device.(L, Cb1, Ct1, 1e11, [3, 1, 1/3] ./ (d0.s1 / d0.L), 1e2)
x0 = 1/3
a = 1/10
for ω in range(1e-3, 6, length=512)
	@printf(file, "%.6g\t%.6g\t%.6g\t%.6g\n", ω,
		ΔU(d[1], ω*d[1].ω1, [x0-a/2, a]),
		ΔU(d[2], ω*d[2].ω1, [x0-a/2, a]),
		ΔU(d[3], ω*d[3].ω1, [x0-a/2, a]))
end
close(file)

file = open("pwmap", "w")
d = semicond_device(L, Cb1, Ct1, 1e11, 50e-12, 1e2)
x0 = 1/2
for ω in range(1e-3, 6, length=512)
	for a in range(1e-3, 1e0, length=512)
		@printf(file, "%.6g\t%.6g\t%.6g\n", ω, a,
			ΔU(d, ω*d.ω1, [x0*(1-a), a]))
	end
	@printf(file, "\n")
end
close(file)

file = open("pwcut", "w")
d = semicond_device(L, Cb1, Ct1, 1e11, 5e-11, 1e2)
for ω in range(1e-3, 8, length=512)
	@printf(file, "%.6g\t%.6g\t%.6g\n", ω,
		ΔU(d, ω*d.ω2, [0, 1]),
		ΔU(d, ω*d.ω2, [.4, .2]))
end
close(file)

file = open("lwmap", "w")
d = graphene_device(L, Cb2, Ct2, 1e12, 5e-12, 5e2)
x0 = 1/2
for ω in range(1e-3, 6, length=512)
	for a in range(1e-3, 1e0, length=512)
		@printf(file, "%.6g\t%.6g\t%.6g\n", ω, a,
			ΔU(d, ω*d.ω1, [x0*(1-a), a]))
	end
	@printf(file, "\n")
end
close(file)

file = open("lwcut", "w")
d = graphene_device(L, Cb2, Ct2, 1e12, 5e-12, 5e2)
for ω in range(1e-3, 8, length=512)
	@printf(file, "%.6g\t%.6g\t%.6g\n", ω,
		ΔU(d, ω*d.ω2, [0, 1]),
		ΔU(d, ω*d.ω2, [.4, .2]))
end
close(file)

file = open("px0map", "w")
d = semicond_device(L, Cb1, Ct1, 1e11, 5e-11, 1e2)
a = 1/10
for ω in range(1e-3, 6, length=512)
	for x0 in range(a/2, 1-a/2, length=512)
		@printf(file, "%.6g\t%.6g\t%.6g\n", ω, x0,
			(ω*d.ω1*a*d.L/d.s1)^-2 * ΔU(d, ω*d.ω1, [x0-a/2, a]))
	end
	@printf(file, "\n")
end
close(file)

file = open("px0x1map", "w")
d = semicond_device(L, Cb1, Ct1, 1e11, 5e-11, 1e2)
x0 = 1/2
a = 1/10
for ω in range(1e-3, 6, length=512)
	for l = range(0, 1-2a, length=256)
		@printf(file, "%.6g\t%.6g\t%.6g\n", ω, l + a,
			ΔU(d, ω*d.ω1, [(x0-a)*(1-l/(1-2a)), a, l, a]))
	end
	@printf(file, "\n")
end
close(file)

