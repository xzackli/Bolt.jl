using Revise
using Bolt
using ForwardDiff
# using PyPlot
using Plots
using BenchmarkTools

par = CosmoParams()
bg = Background(par)
typeof(Bolt.Peebles())
isa(Bolt.Peebles(),IonizationIntegrator)
isa(par,AbstractCosmoParams)
isa(bg,AbstractBackground)
typeof(par)
typeof(bg)
𝕡𝕚=Bolt.PeeblesI(bg,par)
isa(Bolt.RECFAST(bg=bg),IonizationIntegrator)
ih = IonizationHistory(𝕡𝕚, par, bg)

function testih(Ω_b::DT) where DT
    println("omgeab ", Ω_b)
    par = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(par)
    # 𝕡𝕚=Bolt.PeeblesI(bg,par)
    𝕡𝕚=Bolt.Peebles()
    ih = IonizationHistory(𝕡𝕚, par, bg)
    return ih.Xₑ(-5)
end


f(Ω_b) = testih(Ω_b)
#@time
cl = f(0.046)
#@time
∂cl = ForwardDiff.derivative(f, 0.046)

# using Parameters
# @with_kw struct TEST{T} <: IonizationIntegrator @deftype T
#     x::Real
#     a = 1.0
# end
#
# y = TEST(x=1.0)
# function g(y::DT) where DT
#     return sin(y.x)
# end
# g(y)
# ForwardDiff.derivative(g, y)


k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
sf = source_grid(par, bg, ih, k_grid, BasicNewtonian())

ells = 100:50:1200
cl = cltt(ells, par, bg, ih, sf)

# clf()
# plt.plot(ells, cl .* ells.^2, "-")
plot(ells, cl .* ells.^2)
ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
# gcf()
