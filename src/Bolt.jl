module Bolt

export CosmoParams, AbstractCosmoParams
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


abstract type AbstractCosmoParams{T, DT} end

# natural units, stripped eV
@with_kw struct CosmoParams{T} <: AbstractCosmoParams{T, T} @deftype T
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



include("background.jl")
include("recombination.jl")
include("perturbations.jl")

end
