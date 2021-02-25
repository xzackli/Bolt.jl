
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
    Ω_ν = par.Σm_ν*νfac/par.h^2
    return 1 - (par.Ω_r*(1+(7par.N_ν/8)*(4/11)^(4/3))  # dark energy density
                                         + par.Ω_b + par.Ω_m
                                         + Ω_ν) #assume massive nus are non-rel today
end

# Hubble parameter ȧ/a in Friedmann background
function H_a(a, par::AbstractCosmoParams)
    m=par.Σm_ν #FIXME: alllow multiple species
    #use the relation π²/15 Tγ⁴ = ργ = Ωr ρcrit' - BE integral
    Tγ = (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4) #this is the same as 2.725 in eV for default Ω_r
    Tν = (4/11)^(1/3) * Tγ #assume instant decouple for now
    #^This is a bit backwards, should we put TCMB as input?
    # nγ0 = ζ  * 2 / π^2 * Tγ^3 #photon number density today
    # nν0 = 3/11 * nγ0 #neutrino number density (single species) today
    #F = 1.
    #νfac = (90*ζ /(11*π^4)) * (par.Ω_r * par.h^2/ Tγ)#the factor that goes into nr approx to neutrino energy density
    ρ_ν,P_ν = ρP_0(a,par)#m,Tν,F,false)
    ρ_νr = 3P_ν
    ρ_νnr = ρ_ν - ρ_νr
    return H₀(par) * √((par.Ω_m + par.Ω_b ) * a^(-3)
                        + ρ_ν/ρ_crit(par)
                        + par.Ω_r*(1+(7par.N_ν/8)*(4/11)^(4/3)) * a^(-4)
                        + Ω_Λ(par))
end
# conformal time Hubble parameter, aH
ℋ_a(a, par::AbstractCosmoParams) = a * H_a(a, par)

# functions in terms of x
H(x, par::AbstractCosmoParams) = H_a(x2a(x),par)
ℋ(x, par::AbstractCosmoParams) = ℋ_a(x2a(x), par)

# conformal time
function η(x, par::AbstractCosmoParams)
    return quadgk(a -> 1.0 / (a * ℋ_a(a, par)), 0.0, x2a(x),rtol=1e-2)[1] #FIXME rtol
end

#background FD phase space
function f0(q,a,par::AbstractCosmoParams)
    Tν =  (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4) ##assume instant decouple for now
    m = par.Σm_ν  #FIXME allow for multiple species
    T_dec = 1e6 #assume decoupling of 1MeV - this is not quite right but for now
    #Below is fine since if photons/ur are passed m=0, and if massive use neutrio temp as it should
    a_dec= Tν/T_dec #factor of 3 in here somewhere? depends on how many neutrinos have mass?
    (a<a_dec) ? (a_dec = a) : (a_dec= Tν/T_dec) # if we are early enough use the real scale factor
    ϵ = √( q^2 + (a_dec *m)^2) #this is technically right but could just make it q for neutrinos
    gs =  2 #should be 2 for EACH neutrino family (mass eigenstate)
    # if b    #temporary thing to check quadrature with photon integrals
    #     return F* gs / (2π)^3 / ( exp(ϵ/T) -1)
    # else
    return gs / (2π)^3 / ( exp(ϵ/Tν) +1)
    #end

end

function dlnf0dlnq(q,a,par::AbstractCosmoParams) #this is actually only used in perts
    Tν =  (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4) ##assume instant decouple for now
    m = par.Σm_ν  #FIXME allow for multiple species
    T_dec = 1e6 #assume decoupling of 1MeV - this is not quite right but for now
    a_dec= Tν/T_dec
    (a<a_dec) ? (a_dec = a) : (a_dec= Tν/T_dec) # if we are early enough use the real scale factor
    ϵ = √( q^2 + (a_dec * m)^2)
    return q^2 /(ϵ * Tν) /(1 + exp(-ϵ/Tν))
end

#in natural units, q should be 1 at c since c=1 no? will assume this here
#FIXME better quadrature, other codes use asymptotic expansion
function ρP_0(a,par::AbstractCosmoParams)
    #Background phase space energy density and pressure
    Tν =  (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4) ##assume instant decouple for now
    m = par.Σm_ν
    qmin=1e-18 #numerical issue if qmin is smaller - how to choose?
    qmax=1e1 #how to determine qmax?
    #FIXME cheap rtol
    ρ = 4π * a^(-4) * quadgk(q ->  q^2 * √( q^2 + (a*m)^2 ) * f0(q,a,par) ,qmin, qmax,rtol=1e-2)[1]
    P = 4π/3 * a^(-4) * quadgk(q -> q^2 * q^2 /√( q^2 + (a*m)^2) * f0(q,a,par), qmin, qmax,rtol=1e-2)[1]
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
    ℋ::IT
    ℋ′::IT
    ℋ′′::IT
    η::IT
    η′::IT
    η′′::IT
end

function Background(par::AbstractCosmoParams{T}; x_grid=-20.0:0.01:0.0) where T
    ℋ_ = spline(x_grid, [ℋ(x, par) for x in x_grid])
    η_ = spline(x_grid, [η(x, par) for x in x_grid])
    #f0_ = spline()

    return Background(
        T(H₀(par)),
        T(η(0.0, par)),
        T(ρ_crit(par)),
        T(Ω_Λ(par)),
        x_grid,

        ℋ_,
        spline_∂ₓ(ℋ_, x_grid),
        spline_∂ₓ²(ℋ_, x_grid),

        η_,
        spline_∂ₓ(η_, x_grid),
        spline_∂ₓ²(η_, x_grid),
    )
end
