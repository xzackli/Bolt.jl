using Bolt
using ForwardDiff
using PyPlot
using BenchmarkTools

𝕡 = ΛCDMParams()
bg = Background(𝕡)
ih = IonizationHistory(Peebles(), 𝕡, bg)

k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)

sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())



ells = 100:50:1200
cl = cltt(ells, par, bg, ih, sf)

clf()
plt.plot(ells, cl .* ells.^2, "-")
ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
gcf()

error()