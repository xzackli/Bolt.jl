using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using Interpolations

par = CosmoParams()
bg = Background(par)
ih = SahaPeeblesHistory(par, bg)
x_grid = bg.x_grid

∂ₓ(f, x) = Interpolations.gradient(f, x)[1]
∂ₓ²(f, x) = Interpolations.hessian(f, x)[1]

clf()

fig, ax = subplots(1,2,figsize=(10,5))
ax[1].plot(x_grid, ih.τ.(x_grid), "-", label=raw"$\tau$")
ax[1].plot(x_grid, [abs(∂ₓ(ih.τ,x)) for x in x_grid], "--", label=raw"$|\tau^\prime|$")

ax[2].plot(x_grid, [ih.g̃(x) for x in x_grid], "-", label=raw"$\tilde{g}$")
ax[2].plot(x_grid, [∂ₓ(ih.g̃,x) / 10 for x in x_grid], "--", label=raw"$\tilde{g}\prime/10$")
ax[2].plot(x_grid, [∂ₓ²(ih.g̃,x) / 300 for x in x_grid], "--", label=raw"$\tilde{g}\prime/300$")

ax[1].set_yscale("log")
ax[1].legend()
ax[1].set_xlabel(raw"$x$")
ax[2].set_xlim(-8.0, -6.0)
ax[2].set_ylim(-3.5, 5.5)
ax[2].legend()
ax[2].set_xlabel(raw"$x$")
tight_layout()
gcf()
