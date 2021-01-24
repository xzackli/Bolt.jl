using PyPlot
using FFTW
using Bolt
using LinearAlgebra

N = 64
μ = 0
q = 0.0
r = 10 .^ range(-4., 4., length=N)

pl = Bolt.plan_fftlog(r, μ, q, 1.0; kropt=true)
aₙ = r .^ (μ + 1) .* exp.(-r.^2 / 2)
# aₙ = exp.(-(log.(r) ).^2 / 2)
y = similar(r, ComplexF64)
y2 = similar(r, ComplexF64)
mul!(y, pl, aₙ)

clf()
plt.plot(r, aₙ)
# plt.plot(r, aₙ .* r.^(-q))
# ylim(1e-5, 1.2)
plt.plot(pl.k, y)
yscale("symlog", linthreshy=1e-5)
# yscale("log")
xscale("log")
gcf()

##

# y .*= sin.(0.5 .* r)
ldiv!(y2,pl,y)

clf()
plt.plot(r, aₙ)
plt.plot(pl.r, y2)
# yscale("symlog", linthreshy=1e-5)
xscale("log")
gcf()

##
aₙ = sf[1200:end,1]
y = similar(r, ComplexF64)
r = x2a.(bg.x_grid[1200:end])
μ = 30
q = 0.9
pl = Bolt.plan_fftlog(r, μ, q, 1.0; kropt=true)
mul!(y, pl, aₙ)

clf()
plt.plot(r, aₙ)
plt.plot(r, y)

yscale("symlog", linthreshy=maximum(norm.(y)) * 1e-2)
xscale("log")
gcf()


##
using DelimitedFiles
fftdata = readdlm("test/fftlog64.txt", ' ', Float64, '\n')
k_ref = fftdata[:,1]
f_ref = fftdata[:,2]

##
clf()
plot(k_ref .- fftlogfreq(pl))
gcf()
##

isapprox(y, f_ref)

##
clf()
plot(y .- f_ref)
gcf()

##
clf()
plt.plot(fftlogfreq(pl), reverse(y))
plt.plot(k_ref, f_ref)
yscale("symlog", base=10, linthresh=1e-5)
xscale("log")
gcf()
