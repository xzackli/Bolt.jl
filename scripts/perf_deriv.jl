using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using BenchmarkTools
using ForwardDiff
using QuadGK
using Zygote: @adjoint

##
# using Interpolations, Zygote
# itp = interpolate(([0.0,1.0],), [1.0, 2.0], Gridded(Linear()))
# f(c) = Zygote.forwarddiff(c -> itp[c], c)
# # f(c) = itp(c)
# f'(0.5)
using BenchmarkTools
f(x::Vector) = sum(sin, x) + prod(tan, x) * sum(sqrt, x)

##


function test()
    par = CosmoParams()
    xgrid = collect(-18:0.01:0.0)
    zgrid = x2z.(xgrid)
    Xₑ = Bolt.saha_peebles_recombination(par)
    τ = Bolt.τ_function(xgrid, Xₑ, par)
    g̃ = Bolt.g̃_function(par, Xₑ, τ)
    η = Bolt.η_function(xgrid, par)
    # g̃′ = ForwardDiff.gradient()
    # Xₑ'(-1.0)
    # τ'(-1.0) / Bolt.τ′(-1.0, Xₑ, par)

    clf()
    plot(xgrid, [-τ'(x) for x in xgrid], "-")
    plot(xgrid, [-Bolt.τ′(x, Xₑ, par) for x in xgrid], "--")
    # yscale("log")
    # xlim(-10, -9.8)
    # ylim(0.5e4,1.2e4)
    gcf()

    # @btime g̃(x) setup=(x=-10rand())
    # @btime g̃'(x) setup=(x=-10rand())
end

test()
