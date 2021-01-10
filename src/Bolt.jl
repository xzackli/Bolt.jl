module Bolt

export Cosmo

using Parameters
using Unitful, UnitfulAstro, NaturallyUnitful
using NLsolve
using OrdinaryDiffEq

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
    H₀ = ustrip(natural(h * 100.0u"km/s/Mpc"))  # Hubble constant [s^-1]
    ρ_c = (3 / 8π) * H₀^2 /
        ustrip(natural(float(NewtonianConstantOfGravitation))) # critical density [eV⁴]

    # recombination parameters
    Λ_2s_to_1s = ustrip(natural(8.227u"s^-1"))  # rate of hydrogen double transition from 2s to 1s
    ε₀_H = ustrip(natural(13.605698u"eV"))  # ionization energy of hydrogen
    m_e = ustrip(natural(float(ElectronMass)))
    m_H = ustrip(natural(float(ProtonMass)))
    α = ustrip(natural(float(FineStructureConstant)))
end

include("recombination.jl")


end
