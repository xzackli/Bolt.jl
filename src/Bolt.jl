module Bolt

export CosmoParams, AbstractCosmoParams
export Background, AbstractBackground
export IonizationHistory, AbstractIonizationHistory, IonizationIntegrator
export Peebles, PeeblesI
export ρ_σ,ρP_0,f0,dlnf0dlnq,θ,oldH_a #FIXME: quick hack to look at perts
export Hierarchy, boltsolve, BasicNewtonian,unpack,rsa_perts!,boltsolve_rsa
export IE,initial_conditions,unpack,ie_unpack
export source_grid, quadratic_k, cltt,log10_k,plin
export z2a, a2z, x2a, a2x, z2x, x2z, to_ui, from_ui, dxdq

using Parameters
using Unitful, UnitfulAstro
using NLsolve
using OrdinaryDiffEq
using Interpolations
using OffsetArrays
using QuadGK
using ThreadPools
using ForwardDiff, DiffResults
using NumericalIntegration
using FastGaussQuadrature
using StaticArrays
using DoubleFloats
using MuladdMacro
using LinearAlgebra


using FFTW
import SpecialFunctions: logabsgamma, loggamma, gamma, besselj, sphericalbesselj
import HypergeometricFunctions: pochhammer, errcheck
import AbstractFFTs: fftfreq, Plan, plan_fft!, plan_ifft!
import LinearAlgebra: mul!, ldiv!

import UnitfulCosmo, NaturallyUnitful

import PhysicalConstants.CODATA2018: ElectronMass, ProtonMass,
    FineStructureConstant, ThomsonCrossSection, NewtonianConstantOfGravitation


# set the unit system, this should be improved
natural(x) = UnitfulCosmo.mpc(x)
unnatural(x, y) = UnitfulCosmo.unmpc(x, y)
# natural(x) = NaturallyUnitful.natural(x)
# unnatural(x, y) = NaturallyUnitful.unnatural(x, y)

# all unit conversions. should distribute these in-situ someday. Mpc units
const km_s_Mpc_100 = ustrip(natural(100.0u"km/s/Mpc"))  # [Mpc^-1]
const G_natural = ustrip(natural(float(NewtonianConstantOfGravitation))) # [Mpc^2]

abstract type AbstractCosmoParams{T} end

@with_kw struct CosmoParams{T} <: AbstractCosmoParams{T} @deftype T
    h = 0.7  # hubble factor
    Ω_r = 5.0469e-5  # radiation density
    Ω_b = 0.046  # baryon density
    Ω_c = 0.224  # cdm density
    A = 2.097e-9 # scalar amplitude, 1e-10*exp(3.043)
    n = 1.0  # scalar spectral index
    Y_p = 0.24  # primordial helium fraction
    N_ν = 3.046 #effective number of relativisic species (PDG25 value)
    Σm_ν = 0.06 #sum of neutrino masses (eV), Planck 15 default ΛCDM value
end

include("util.jl")
include("bessel/weniger.jl")
include("bessel/moments.jl")
include("background.jl")
include("ionization/ionization.jl")
include("ionization/recfast.jl")
include("perturbations.jl")
include("spectra.jl")

end
