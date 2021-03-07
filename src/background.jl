
# NOTE: Bolt's background functions are in terms of x ≡ ln(a), the log scale factor

# derived quantities (I've chosen natural units, possibly the wrong choice)
const km_s_Mpc_100 = ustrip(natural(100.0u"km/s/Mpc"))  # [eV]
const G_natural = ustrip(natural(float(NewtonianConstantOfGravitation))) #[eV^-2]
const ζ = 1.2020569 #Riemann ζ(3) for phase space integrals
#In the natural units used here, ħ=kb=c=1, G as above

H₀(par::AbstractCosmoParams) = par.h * km_s_Mpc_100
ρ_crit(par::AbstractCosmoParams) = (3 / 8π) * H₀(par)^2 / G_natural  # [eV⁴]
function Ω_Λ(par::AbstractCosmoParams)
    #Below can definitely be more streamlined, I am just making it work for now
    Tγ = (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    νfac = (90 * ζ /(11 * π^4)) * (par.Ω_r * par.h^2 / Tγ)#the factor that goes into nr approx to neutrino energy density
    Ω_ν = par.Σm_ν*νfac/par.h^2 #FIXME I think this is right for a single neutrino? May be mssing Neff/3 factor
    return 1 - (par.Ω_r*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3))  # dark energy density
                                         + par.Ω_b + par.Ω_m
                                         + Ω_ν) #assume massive nus are non-rel today
end

# Hubble parameter ȧ/a in Friedmann background
function H_a(a, par::AbstractCosmoParams)
    ρ_ν,_ = ρP_0(a,par) # we don't atually need pressure?
    return H₀(par) * √((par.Ω_m + par.Ω_b ) * a^(-3)
                        + ρ_ν/ρ_crit(par)
                        + par.Ω_r*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3)) * a^(-4)
                        + Ω_Λ(par))
end
# conformal time Hubble parameter, aH
ℋ_a(a, par::AbstractCosmoParams) = a * H_a(a, par)

# functions in terms of x
H(x, par::AbstractCosmoParams) = H_a(x2a(x),par)
ℋ(x, par::AbstractCosmoParams) = ℋ_a(x2a(x), par)

# conformal time
function η(x, par::AbstractCosmoParams)
    return quadgk(a -> 1.0 / (a * ℋ_a(a, par)), 0.0, x2a(x),rtol=1e-6)[1]
end

#background FD phase space
function f0(q,par::AbstractCosmoParams)
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    #m = par.Σm_ν  #FIXME allow for multiple species
    gs =  2 #should be 2 for EACH neutrino family (mass eigenstate)
    return gs / (2π)^3 / ( exp(q/Tν) +1)
end

function dlnf0dlnq(q,par::AbstractCosmoParams) #this is actually only used in perts
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4) ##assume instant decouple for now
    #m = par.Σm_ν  #FIXME allow for multiple species
    return -q / Tν /(1 + exp(-q/Tν))
end

#in natural units, q should be 1 at c since c=1 no? will assume this here
#FIXME better quadrature, other codes use asymptotic expansion
function ρP_0(a,par::AbstractCosmoParams)
    #Background phase space energy density and pressure
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4) ##assume instant decouple for now
    m = par.Σm_ν
    qmin=1e-18 #numerical issue if qmin is smaller - how to choose?
    qmax=1e1 #how to determine qmax?
    #FIXME cheap rtol
    ρ = 4π * a^(-4) * quadgk(q ->  q^2 * √( q^2 + (a*m)^2 ) * f0(q,par) ,qmin, qmax,rtol=1e-6)[1]
    P = 4π/3 * a^(-4) * quadgk(q -> q^2 * q^2 /√( q^2 + (a*m)^2) * f0(q,par), qmin, qmax,rtol=1e-6)[1]
    return ρ,P#,norm
end
# now build a Background with these functions

# a background is parametrized on the scalar type T, the interpolator type IT,
# and a type for the grid GT
abstract type AbstractBackground{T, IT<:AbstractInterpolation{T,1}, GT} end

struct Background{T, IT, GT} <: AbstractBackground{T, IT, GT}
    H₀::T
    η₀::T
    ρ_crit::T
    Ω_Λ::T

    x_grid::GT
#    logq_grid::GT
    quad_pts::Array{T,1}
    quad_wts::Array{T,1}

    ℋ::IT
    ℋ′::IT
    ℋ′′::IT
    η::IT
    η′::IT
    η′′::IT


    # f0::IT
    # df0::IT
end

function Background(par::AbstractCosmoParams{T}; x_grid=-20.0:0.01:0.0, nq=15) where T
                    #,logq_grid=-6.0:0.1:1.0) where T
    ℋ_  = spline(x_grid, [ℋ(x, par) for x in x_grid])
    η_   = spline(x_grid, [η(x, par) for x in x_grid])
    quad_pts, quad_wts =  gausslegendre( nq ) #12 should get 1e-3, 15 conservative for 1e-6
    #logq_grid   #probably bad but just to get started
    # f0_  = spline(logq_grid, [f0(10.0 ^(lq),par) for lq in logq_grid])
    # df0_ = spline(logq_grid, [dlnf0dlnq(10.0 ^(lq),par) for lq in logq_grid]) #maybe better to do spline gradient?
    return Background(
        T(H₀(par)),
        T(η(0.0, par)),
        T(ρ_crit(par)),
        T(Ω_Λ(par)),

        x_grid,
        #logq_grid,
        convert(Array{T,1},quad_pts), #explicit call to convert instead of constructor for arrays
        convert(Array{T,1},quad_wts),

        ℋ_,
        spline_∂ₓ(ℋ_, x_grid),
        spline_∂ₓ²(ℋ_, x_grid),

        η_,
        spline_∂ₓ(η_, x_grid),
        spline_∂ₓ²(η_, x_grid),


        # f0_,
        # df0_,
    )
end
