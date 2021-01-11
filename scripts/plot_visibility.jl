using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful

par = Cosmo()
xgrid = collect(-18:0.02:0.0)
zgrid = x2z.(xgrid)
Xₑ = Bolt.saha_peebles_recombination(par)
τ = Bolt.τ_function(xgrid, Xₑ, par)
g̃ = Bolt.g̃_function(par, Xₑ, τ)

clf()
fig, axes = subplots(1, 2,figsize=(10,5))
axes[1].plot(xgrid, τ.(xgrid), "-", label=raw"$\tau$")
axes[1].plot(xgrid, [-1 * Bolt.τ′(x, Xₑ, par) for x in xgrid], "--", label=raw"$|\tau^\prime|$")
axes[1].set_yscale("log")
axes[1].legend()
axes[1].set_xlabel(raw"$x$")

axes[2].plot(xgrid, g̃.(xgrid), "-", label=raw"$\tilde{g}$")
axes[2].set_xlim(-8.0, -6.0)
axes[2].legend()
axes[2].set_xlabel(raw"$x$")
tight_layout()
gcf()
