module Bolt

export Cosmo

using Parameters
using Unitful, UnitfulAstro, NaturallyUnitful
using NLsolve
using DifferentialEquations

import PhysicalConstants.CODATA2018: ElectronMass, ProtonMass,
    FineStructureConstant, ThomsonCrossSection, NewtonianConstantOfGravitation

abstract type AbstractCosmo{T} end

# natural units, stripped eV
@with_kw struct Cosmo{T} <: AbstractCosmo{T} @deftype T

    # cosmological parameters
    h = 0.7
    Ω_r = 5.042e-5  # radiation density
    Ω_b = 0.046  # baryon density
    Ω_m = 0.224  # matter density
    Ω_Λ = 1 - (Ω_r + Ω_b + Ω_m)  # dark energy density
    n = 1.0  # spectral index
    Y_p = 0.0  # primordial helium fraction
    T₀ = ustrip(natural(2.725u"K"))  # CMB temperature [K]

    # derived parameters
    H₀ = ustrip(natural(h * 100.0u"km/s/Mpc"))  # Hubble constant [s^-1]
    ρ_c = (3 / 8π) * H₀^2 /
        ustrip(natural(float(NewtonianConstantOfGravitation))) # critical density [eV⁴]

end

# Hubble parameter ȧ/a in Friedmann background
H(a, par) = par.H₀ * √((par.Ω_m + par.Ω_b) * a^(-3) + par.Ω_r * a^(-4) + par.Ω_Λ)

# conformal time Hubble parameter, (1/a) * da/dη
ℋ(a, par) = par.H₀ * √((par.Ω_m + par.Ω_b) * a^(-1) + par.Ω_r * a^(-2) + par.Ω_Λ * a^2)

# note: x ≡ ln a



include("recombination.jl")


end
