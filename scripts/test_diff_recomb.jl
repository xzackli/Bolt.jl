# test that the ionization history is differentiable

using Bolt
using Parameters
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using ForwardDiff

@with_kw struct GradientCosmo{T, DT} <: Bolt.AbstractCosmo{T, DT} @deftype T
    h = 0.7
    Ω_r = 5.042e-5  # radiation density
    Ω_b::DT = 0.046  # baryon density
    Ω_m = 0.224  # matter density
    n = 1.0  # spectral index
    Y_p = 0.0  # primordial helium fraction
    T₀ = ustrip(natural(2.725u"K"))  # CMB temperature [K]
end

##
function f(Ω_b)
    par = GradientCosmo(Ω_b=Ω_b)
    z_transition = 1587.4
    x_transition = z2x(z_transition)

    Xₑ = Bolt.saha_peebles_recombination(par)
    # return Xₑ(a_transition * 1.2)
    xgrid = collect(0.0:-0.05:-18)
    return (Bolt.τ_x(xgrid, Xₑ, par))[end]
end


ForwardDiff.derivative(f, 0.046)
