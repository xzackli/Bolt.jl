using Bolt
using ForwardDiff
using PyPlot
using BenchmarkTools

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)

k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
@time sf = source_function(par, bg, ih, k_grid, BasicNewtonian())

##
clf()
k̂ = 340bg.H₀
bes = Bolt.bessel_interpolator(100, k_grid[end] * bg.η₀)
plot(bg.x_grid, [sf(x, k̂) * bes(k̂*abs(bg.η₀ - bg.η(x)))/1e-3 for x in bg.x_grid], "-", lw=0.5)

ylabel(raw"Source function $\times$ bessel / $10^{-3}$")
xlabel(raw"$x$")
xlim(-8, 0)
ylim(-1, 3.5)
# xlim(-8,-6)
gcf()
##


dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000)
cltt(1000, sf, dense_kgrid, par, bg)

##

using ThreadPools
ThreadPools.qmap()

ells = 100:20:1200
@time cl = thCl(ells,  sf, dense_kgrid, par, bg)

##
clf()
plt.plot(ells, cl .* ells.^2, "-")
ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
# yscale("log")
gcf()

##
