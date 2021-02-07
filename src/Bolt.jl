module Bolt

export CosmoParams, AbstractCosmoParams
export Background, AbstractBackground
export IonizationHistory, AbstractIonizationHistory
export Peebles
export Hierarchy, boltsolve, BasicNewtonian
export source_grid, quadratic_k, cltt
export z2a, a2z, x2a, a2x, z2x, x2z

using Parameters
using Unitful, UnitfulAstro, NaturallyUnitful
using NLsolve
using OrdinaryDiffEq
using Interpolations
using OffsetArrays
using QuadGK
using ThreadPools
using ForwardDiff, DiffResults
using NumericalIntegration

using FFTW
import SpecialFunctions: lgamma, sphericalbesselj
import AbstractFFTs: fftfreq, Plan, plan_fft!, plan_ifft!
import LinearAlgebra: mul!, ldiv!


import PhysicalConstants.CODATA2018: ElectronMass, ProtonMass,
    FineStructureConstant, ThomsonCrossSection, NewtonianConstantOfGravitation


abstract type AbstractCosmoParams{T} end

@with_kw struct CosmoParams{T} <: AbstractCosmoParams{T} @deftype T
    h = 0.7  # hubble factor
    Ω_r = 5.042e-5  # radiation density
    Ω_b = 0.046  # baryon density
    Ω_m = 0.224  # matter density
    n = 1.0  # spectral index
    Y_p = 0.0  # primordial helium fraction, currently unused
end

include("util.jl")
include("fftlog.jl")
include("background.jl")
include("ionization.jl")
include("perturbations.jl")
include("spectra.jl")

end
