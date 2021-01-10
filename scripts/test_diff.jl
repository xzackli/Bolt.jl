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
    a_transition = z2a(z_transition)
    saha_z_grid = 1800:-10:z_transition
    peebles_z_grid = z_transition:-10:100
    early_time_Xₑ = Bolt.saha_Xₑ(par)
    return early_time_Xₑ(a_transition)
end


ForwardDiff.derivative(f, 0.046)
