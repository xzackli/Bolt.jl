module Bolt

export Cosmo, AbstractCosmo
export z2a, a2z, x2a, a2x, z2x, x2z

using Parameters
using Unitful, UnitfulAstro, NaturallyUnitful
using NLsolve
using OrdinaryDiffEq
using NumericalIntegration
using Interpolations
using ForwardDiff
using OffsetArrays
using QuadGK
using SpecialFunctions
using ThreadPools
# using Zygote

import PhysicalConstants.CODATA2018: ElectronMass, ProtonMass,
    FineStructureConstant, ThomsonCrossSection, NewtonianConstantOfGravitation

abstract type AbstractCosmo{T, DT} end

# natural units, stripped eV
@with_kw struct Cosmo{T} <: AbstractCosmo{T, T} @deftype T
    h = 0.7  # hubble factor
    Ω_r = 5.042e-5  # radiation density
    Ω_b = 0.046  # baryon density
    Ω_m = 0.224  # matter density
    n = 1.0  # spectral index
    Y_p = 0.0  # primordial helium fraction
    T₀ = ustrip(natural(2.725u"K"))  # CMB temperature [K]

    lmax::Int = 3
    ℓᵧ::Int = 30  # Boltzmann hierarchy cutoff
end

# # utility function to make scalar interpolators
# function differentiable_interpolator(xgrid, ygrid)
#     itp = interpolate((xgrid,), ygrid, Gridded(Linear()))
#     f(c) = Zygote.forwarddiff(c -> itp(c), c)
#     return f
# end

# derived quantities (make sure to convert to natural units)
const km_s_Mpc_100 = ustrip(natural(100.0u"km/s/Mpc"))  # [eV]
H₀(par::AbstractCosmo) = par.h * km_s_Mpc_100
const G_natural = ustrip(natural(float(NewtonianConstantOfGravitation)))
ρ_crit(par::AbstractCosmo) = (3 / 8π) * H₀(par)^2 / G_natural  # [eV⁴]
Ω_Λ(par::AbstractCosmo) = 1 - (par.Ω_r + par.Ω_b + par.Ω_m)  # dark energy density

# Hubble parameter ȧ/a in Friedmann background
H_a(a, par) = H₀(par) * √((par.Ω_m + par.Ω_b) * a^(-3) + par.Ω_r * a^(-4) + Ω_Λ(par))
H(x, par) = H_a(x2a(x),par)

# conformal time Hubble parameter, aH
ℋ_a(a, par) = a * H_a(a, par)
ℋ(x, par) = ℋ_a(x2a(x), par)

function ℋ′(x, par)
    a = x2a(x)
    return -H₀(par) * (2par.Ω_r + (par.Ω_b + par.Ω_m) * a - 2Ω_Λ(par) * a^4) /
        (2 * a * √(par.Ω_r + (par.Ω_b + par.Ω_m) * a + Ω_Λ(par) * a^4))
end


# conformal time
function η_function(xgrid, par::AbstractCosmo{T,DT}) where {T, DT}
    agrid = x2a.(xgrid)
    η = [quadgk(ap -> 1.0 / (ap * Bolt.ℋ_a(ap, par)), 0.0, a)[1] for a in agrid]
    return interpolate((xgrid,), η, Gridded(Linear()))
end



# utilities for x ↔ scale factor ↔ redshift
a2z(a::T) where T = one(T)/a - one(T)
z2a(z::T) where T = one(T)/(z + one(T))
a2x(a) = log(a)
x2a(x) = exp(x)
z2x(z) = a2x(z2a(z))
x2z(x) = a2z(x2a(x))

include("recombination.jl")
include("perturbations.jl")

end
