module Bolt

export CosmoParams, AbstractCosmoParams
export Background, AbstractBackground
export IonizationHistory, AbstractIonizationHistory, IonizationIntegrator
export Peebles, PeeblesI
export ρ_σ,ρP_0,f0,dlnf0dlnq,θ,oldH_a #FIXME: quick hack to look at perts
export Hierarchy, boltsolve, BasicNewtonian,unpack,rsa_perts!,boltsolve_rsa
export source_grid, quadratic_k, cltt,log10_k,plin
export z2a, a2z, x2a, a2x, z2x, x2z, to_ui, from_ui, dxdq

using ComponentArrays
using FastGaussQuadrature
using ForwardDiff, DiffResults
using Interpolations
using NLsolve
using NumericalIntegration
using OffsetArrays: OffsetArray, Origin
using OrdinaryDiffEq
using Parameters
using QuadGK
using ThreadPools
using Unitful, UnitfulAstro, NaturallyUnitful

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
    A = 2.097e-9 # scalar amplitude, 1e-10*exp(3.043)
    n = 1.0  # spectral index
    Y_p = 0.24  # primordial helium fraction
    N_ν = 3.046 #effective number of relativisic species (PDG25 value)
    Σm_ν = 0.06 #sum of neutrino masses (eV), Planck 15 default ΛCDM value
end

include("util.jl")
include("background.jl")
include("ionization/ionization.jl")
include("ionization/recfast.jl")
include("perturbations.jl")
include("spectra.jl")


end
