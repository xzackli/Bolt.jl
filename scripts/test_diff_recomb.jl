# test that the ionization history is differentiable

using Bolt
using Parameters
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using ForwardDiff

par = CosmoParams()
xgrid = collect(-18:0.001:-0.0)
zgrid = x2z.(xgrid)
Xₑ = Bolt.saha_peebles_recombination(par)
τ, τ′, τ′′ = Bolt.τ_functions(xgrid, Xₑ, par)
g̃, g̃′ = Bolt.g̃_functions(τ, τ′, τ′′)

clf()
fig, axes = subplots(1, 2,figsize=(10,5))
axes[1].plot(xgrid, τ.(xgrid), "-", label=raw"$\tau$")
axes[1].plot(xgrid, [-1 * τ′(x) for x in xgrid], "--", label=raw"$|\tau^\prime|$")
axes[1].set_yscale("log")
axes[1].legend()
axes[1].set_xlabel(raw"$x$")

using ForwardDiff
∂ₓ(f, x) = ForwardDiff.derivative(f, x)
g̃′′ = x -> ∂ₓ(g̃′, x)

axes[2].plot(xgrid, g̃.(xgrid), "-", label=raw"$\tilde{g}$")
axes[2].plot(xgrid, g̃′.(xgrid) / 10, "--", label=raw"$\tilde{g}\prime$")
axes[2].plot(xgrid, g̃′′.(xgrid) / 300, ":", label=raw"$\tilde{g}\prime$")

axes[2].set_xlim(-8.0, -6.0)
axes[2].set_ylim(-3.5, 5.5)
axes[2].legend()
axes[2].set_xlabel(raw"$x$")
tight_layout()
gcf()

##
# ∂ₓ²(f, x::T) where T <: Real = ForwardDiff.hessian(x -> f(x[1]), [x])[1]
∂ₓ²(f,x) = ∂ₓ(x->∂ₓ(f,x),x)
# numdiff(f,x) = (f(x+1e-3) - f(x-1e-3)) / 2e-3

clf()
using Zygote
ddx = [∂ₓ(g̃,x) for x in xgrid]

h = x -> ∂ₓ(g̃,x)

plot(xgrid[1:end], ddx, "-", label=raw"$\tilde{g}$")
# plot(xgrid, [h'(x) for x in xgrid], "-", label=raw"$\tilde{g}$")
plot(xgrid, [∂ₓ(g̃,x) for x in xgrid], "--", label=raw"$\tilde{g}$")
xlim(-8,-6)
gcf()
