using Bolt
using ForwardDiff
using PyPlot
using BenchmarkTools

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
sf = sourcefunction(par, bg, ih, k_grid, BasicNewtonian())

using ThreadPools
ells = 100:20:1200
cl = qmap(ℓ->cltt(ℓ, par, bg, ih, sf), ells)
clf()
plt.plot(ells, cl .* ells.^2, "-")
ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
gcf()
