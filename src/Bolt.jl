module Bolt

export CosmoParams, AbstractCosmoParams
export Background, AbstractBackground
export IonizationHistory, AbstractIonizationHistory, IonizationIntegrator
export Peebles, PeeblesI
export ρ_σ,ρP_0,f0,dlnf0dlnq,θ,oldH_a #FIXME: quick hack to look at perts
export Hierarchy, boltsolve, BasicNewtonian,unpack,rsa_perts!,boltsolve_rsa
export source_grid, quadratic_k, cltt,log10_k,plin, source_grid_P, clte,clee
export z2a, a2z, x2a, a2x, z2x, x2z, to_ui, from_ui, dxdq

using Parameters
using Unitful, UnitfulAstro
using NonlinearSolve
using OrdinaryDiffEq
using DataInterpolations
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
import Interpolations
import SpecialFunctions: logabsgamma, loggamma, gamma, sphericalbesselj
import HypergeometricFunctions: pochhammer, errcheck, pFqmaclaurin
import AbstractFFTs: fftfreq, Plan, plan_fft!, plan_ifft!
import LinearAlgebra: mul!, ldiv!
import Bessels: besselj

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
const mass_natural = ustrip(UnitfulCosmo.mpc(1.0u"eV")) # [Mpc^-1] #FIXME? Might want to put neutrino conversion elsewhere, but need to convert eV


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
    Σm_ν = 0.06*mass_natural #sum of neutrino masses (eV), Planck 15 default ΛCDM value
end

include("util.jl")

include("bessel/weniger.jl")
include("bessel/moments.jl")
include("bessel/interpolator.jl")
include("bessel/integrator.jl")

include("background.jl")
include("ionization/ionization.jl")
include("ionization/recfast.jl")
include("perturbations.jl")
include("spectra.jl")

end
