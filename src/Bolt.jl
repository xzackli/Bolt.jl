module Bolt

export Cosmo, AbstractCosmo
export z2a, a2z, x2a, a2x, z2x, x2z

using Parameters
using Unitful, UnitfulAstro, NaturallyUnitful
using NLsolve
using OrdinaryDiffEq
using NumericalIntegration
using Interpolations

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
end

# derived quantities (make sure to convert to natural units)
const km_s_Mpc_100 = ustrip(natural(100.0u"km/s/Mpc"))  # [eV]
H₀(par::AbstractCosmo) = par.h * km_s_Mpc_100
const G_natural = ustrip(natural(float(NewtonianConstantOfGravitation)))
ρ_crit(par::AbstractCosmo) = (3 / 8π) * H₀(par)^2 / G_natural  # [eV⁴]
Ω_Λ(par::AbstractCosmo) = 1 - (par.Ω_r + par.Ω_b + par.Ω_m)  # dark energy density

# Hubble parameter ȧ/a in Friedmann background
H(a, par) = H₀(par) * √((par.Ω_m + par.Ω_b) * a^(-3) + par.Ω_r * a^(-4) + Ω_Λ(par))

# conformal time Hubble parameter, (1/a) * da/dη
ℋ(a, par) = H₀(par) * √((par.Ω_m + par.Ω_b) * a^(-1) + par.Ω_r * a^(-2) + Ω_Λ(par) * a^2)

# conformal time
η(a::T, par) where T = quadgk(a′ -> one(T) / (a′ * ℋ(a′, par)), zero(T), a)[1]

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
