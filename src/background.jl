
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
    #^This is a bit backwards, should we put TCMB as input?
    nγ0 = ζ  * 2 / π^2 * Tγ^3 #photon number density today
    nν0 = 3/11 * nγ0 #neutrino number density (single species) today
    Tν = (4/11)^(1/3) * Tγ #assume instant decouple for now
    F = 1.
    νfac = (90*ζ /(11*π^4)) * (par.Ω_r * par.h^2/ Tγ)#the factor that goes into nr approx to neutrino energy density
    #a=1
    #Find norm so we get the correct number density of neutrinos today
    #ρ_ν,P_ν,norm = ρP_0(a,m,Tν,F,false)
    ρ_ν,P_ν = ρP_0(a,m,Tν,F,false)
    # println("m *nν0: ", m*nν0)# * (8*π /3 *G_natural/km_s_Mpc_100^2) )# (0.00001014051))
    # println("normalized density: ",ρ_ν )
    # println("expected rho: ", 7/8 * 2 * π^2 /30 * Tν^4 )

    #test with photons
    # ργ,Pγ,_ = ρP_0(a,0,Tγ,F,true)
    # ρur,Pur,_ = ρP_0(a,0,Tν,par.N_ν,false) #need neff factor for FD
    ρ_νr = 3P_ν
    ρ_νnr = ρ_ν - ρ_νr
    # println("1/nufac: ",1/νfac)
    # println("mass*nufac: ", m*νfac)
    # println(" a= ",a)
    # println("TCMB: ", Tγ)
    # println("Mnu: ",m)
    # println(" nr piece  ",ρ_νnr, " r piece: ",ρ_νr)
    # println("total neutrino density/ρcrit: ", ρ_ν/ρ_crit(par))
    # println("nr neutrino density/ρcrit: ", ρ_νnr/ρ_crit(par))
    # println(" nr approx to neutrino Ω: ", par.Σm_ν*νfac/par.h^2 * a^(-3))
    # println(" analytic ρur/ρcrit = ", π^2 /15 * (7/8)*par.N_ν*Tν^(4) *a^(-4)  /ρ_crit(par) )
    # println(" analytic ργ/ρcrit = ", π^2 /15 * Tγ^(4) *a^(-4)  /ρ_crit(par) )
    # println(" ργ/(3Pγ) = ", ργ /3. /Pγ) #this should be 1, and it is
    # println(" ργ/ρcrit = ", ργ/ρ_crit(par), " Omegaγ= ", par.Ω_r*a^(-4) ) #but this is not right...
    # println(" ρur/ρcrit = ", ρur/ρ_crit(par), " Omegaur= ", par.Ω_r*((7par.N_ν/8)*(4/11)^(4/3))*a^(-4) ) #but this is not right...
    # println(" ρur/(3Pur) = ", ρur /3. /Pur) #this should be 1,
    # println(" ρν/(3Pν) = ", ρ_ν /3. /P_ν)
    #^the energy density units must be wrong, need to convert into same units as critical density?
    #since it is 1, I would guess problem is 1. unit conversion or 2. f0 normalization (which also may be units)

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
function f0(q,a,m,T,F,b)
    T_dec = 1e6 #assume decoupling of 1MeV - this is not quite right but for now
    #Below is fine since if photons/ur are passed m=0, and if massive use neutrio temp as it should
    a_dec= T/T_dec #factor of 3 in here somewhere? depends on how many neutrinos have mass?
    (a<a_dec) ? (a_dec = a) : (a_dec= T/T_dec) # if we are early enough use the real scale factor
    ϵ = √( q^2 + (a_dec *m)^2) #this is technically right but could just make it q for neutrinos
    #FIXME allow for multiple species
    gs =  2 #should be 2 for EACH neutrino family (mass eigenstate)
    if b    #temporary thing to check quadrature with photon integrals
        return F* gs / (2π)^3 / ( exp(ϵ/T) -1)
    else
        return F* gs / (2π)^3 / ( exp(ϵ/T) +1)
    end

end

function dlnf0dlnq(q,a,m,T) #this is actually only used in perts
    T_dec = 1e6 #assume decoupling of 1MeV - this is not quite right but for now
    a_dec= T/T_dec
    (a<a_dec) ? (a_dec = a) : (a_dec= T/T_dec) # if we are early enough use the real scale factor
    ϵ = √( q^2 + (a_dec * m)^2)
    return q^2 /(ϵ * T) /(1 + exp(-ϵ/T))
end

#in natural units, q should be 1 at c since c=1 no? will assume this here
#FIXME better quadrature, other codes use asymptotic expansion with the spline
function ρP_0(a,m,T,F,b)
    #Background phase space energy density and pressure
    qmin=1e-18 #numerical issue if qmin is smaller - how to choose?
    qmax=1e1 #how to determine qmax?
    #norm = 4π * a^(-3) *quadgk(q ->  f0(q,a,m,T,F,b) ,qmin, qmax)[1]
    #FIXME cheap rtol
    ρ = 4π * a^(-4) * quadgk(q ->  q^2 * √( q^2 + (a*m)^2 ) * f0(q,a,m,T,F,b) ,qmin, qmax,rtol=1e-2)[1]
    P = 4π/3 * a^(-4) * quadgk(q -> q^2 * q^2 /√( q^2 + (a*m)^2) * f0(q,a,m,T,F,b), qmin, qmax,rtol=1e-2)[1]
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
