using PyPlot
using FFTW
using Bolt
using LinearAlgebra

N = 64
μ = 0
q = 0.5
r = 10 .^ range(-4., 4., length=N)

pl = Bolt.plan_fftlog(r, μ, q, 1.0; kropt=true)
aₙ = r .^ (μ + 1) .* exp.(-r.^2 / 2)
y = similar(r, ComplexF64)
y2 = similar(r, ComplexF64)
mul!(y, pl, aₙ)

clf()
plt.plot(r, aₙ)
plt.plot(pl.k, y)
yscale("symlog", linthreshy=1e-5)
xscale("log")
gcf()

##
ldiv!(y2,pl,y)

clf()
plt.plot(r, aₙ)
plt.plot(pl.r, y2)
yscale("symlog", linthreshy=1e-5)
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
