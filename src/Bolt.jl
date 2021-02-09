module Bolt

using AbstractFFTs: fftfreq, Plan, plan_fft!, plan_ifft!
using FFTW
using ForwardDiff, DiffResults
using Interpolations
using NaturallyUnitful
using NLsolve
using NumericalIntegration
using OffsetArrays
using OrdinaryDiffEq
using Parameters
using PhysicalConstants.CODATA2018: ElectronMass, ProtonMass,
    FineStructureConstant, ThomsonCrossSection, NewtonianConstantOfGravitation
using QuadGK
using SpecialFunctions: lgamma, sphericalbesselj
using ThreadPools
using Unitful
using UnitfulAstro

import LinearAlgebra: mul!, ldiv!

export ΛCDMParams, AbstractParams,
    Background, AbstractBackground,
    IonizationHistory, AbstractIonizationHistory,
    Peebles,
    Hierarchy, boltsolve, BasicNewtonian,
    source_grid, quadratic_k, cltt,
    z2a, a2z, x2a, a2x, z2x, x2z


abstract type AbstractParams{T} end

@with_kw struct ΛCDMParams{T} <: AbstractParams{T} @deftype T
    h  = 0.7       # hubble factor
    Ωr = 5.042e-5  # radiation density
    Ωb = 0.046     # baryon density
    Ωm = 0.224     # matter density
    n  = 1.0       # spectral index
    Yp = 0.0       # primordial helium fraction, currently unused
    Nν = 3.046     # effective number of relativisic species (PDG25 value)
end

include("util.jl")
include("constants.jl")
include("fftlog.jl")
include("background.jl")
include("ionization.jl")
include("perturbations.jl")
include("spectra.jl")

end
