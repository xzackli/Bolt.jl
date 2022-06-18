
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
    νfac =(90 * ζ /(11 * π^4)) * (par.Ω_r * par.h^2 / Tγ) *((par.N_ν/3)^(3/4))
    #^the factor that goes into nr approx to neutrino energy density, plus equal sharing ΔN_eff factor for single massive neutrino
    Ω_ν = par.Σm_ν*νfac/par.h^2
    return 1 - (par.Ω_r*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3))  # dark energy density
                                         + par.Ω_b + par.Ω_m
                                         + Ω_ν
                                         ) #assume massive nus are non-rel today
end

#background FD phase space
function f0(q,par::AbstractCosmoParams)
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    gs =  2 #should be 2 for EACH neutrino family (mass eigenstate)
    return gs / (2π)^3 / ( exp(q/Tν) +1)
end

function dlnf0dlnq(q,par::AbstractCosmoParams) #this is actually only used in perts
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    return -q / Tν /(1 + exp(-q/Tν))
end

#This is just copied from perturbations.jl for now - but take out Pressure - maybe later restore for FD tests?
function ρP_0(a,par::AbstractCosmoParams,quad_pts,quad_wts)
    #Do q integrals to get the massive neutrino metric perturbations
    #MB eqn (55)
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    #Not allowed to set Neff=0 o.w. breaks this #FIXME add an error message
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    #FIXME: avoid repeating code? and maybe put general integrals in utils?
    m = par.Σm_ν
    ϵx(x, am) = √(xq2q(x,logqmin,logqmax)^2 + (am)^2)
    Iρ(x) = xq2q(x,logqmin,logqmax)^2  * ϵx(x, a*m) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    IP(x) = xq2q(x,logqmin,logqmax)^2  * (xq2q(x,logqmin,logqmax)^2 /ϵx(x, a*m)) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    xq,wq =quad_pts,quad_wts
    ρ = 4π * a^(-4) * sum(Iρ.(xq).*wq)
    P = 4π/3 * a^(-4) *sum(IP.(xq).*wq)
    return ρ,P
end

#neglect neutrinos, this is for ionization debugging purposes only
function oldH_a(a, par::AbstractCosmoParams)
    return H₀(par) * √((par.Ω_m + par.Ω_b ) * a^(-3)
                        + par.Ω_r*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3)) * a^(-4)
                        + Ω_Λ(par))
end

# Hubble parameter ȧ/a in Friedmann background
function H_a(a, par::AbstractCosmoParams,quad_pts,quad_wts)
    #ρ_ν,_ = ρP_0(a,par,quad_pts,quad_wts) # we don't atually need pressure?
    ρ_ν,_ = ρP_0(a,par,quad_pts,quad_wts) #FIXME dropped pressure, need to decide if we want it for tests?
    #ρ_ν = ρℳ(a2x(a))
    return H₀(par) * √((par.Ω_m + par.Ω_b ) * a^(-3)
                        + ρ_ν/ρ_crit(par)
                        + par.Ω_r* a^(-4)*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3))
                        + Ω_Λ(par))
end
# conformal time Hubble parameter, aH
ℋ_a(a, par::AbstractCosmoParams,quad_pts,quad_wts) = a * H_a(a, par,quad_pts,quad_wts)

# functions in terms of x
H(x, par::AbstractCosmoParams,quad_pts,quad_wts) = H_a(x2a(x),par,quad_pts,quad_wts)
ℋ(x, par::AbstractCosmoParams,quad_pts,quad_wts) = ℋ_a(x2a(x), par,quad_pts,quad_wts)

# conformal time
function η(x, par::AbstractCosmoParams,quad_pts,quad_wts)
    #fast copy from q - need to check accuracy (#FIXME) but a plays the role of q
    logamin,logamax=-11.75,log10(x2a(x)) #0,x2a(x)
    #convert ui to a,for now pick
    Iη(y) = 1.0 / (xq2q(y,logamin,logamax) * ℋ_a(xq2q(y,logamin,logamax), par,quad_pts,quad_wts))/ dxdq(xq2q(y,logamin,logamax),logamin,logamax)
    #return quadgk(a -> 1.0 / (a * ℋ_a(a, par)), 0.0, x2a(x),rtol=1e-6)[1]
    return sum(Iη.(quad_pts).*quad_wts)
end

# now build a Background with these functions

# a background is parametrized on the scalar type T, the interpolator type IT,
# and a type for the grid GT
abstract type AbstractBackground{T, IT<:AbstractInterpolation{T,1}, GT, U} end

struct Background{T, IT, GT, U, AT<:AbstractArray{U,1}} <: AbstractBackground{T, IT, GT, U}
    H₀::T
    η₀::T
    ρ_crit::T
    Ω_Λ::T

    x_grid::GT
    quad_pts::AT
    quad_wts::AT

    ℋ::IT
    ℋ′::IT
    ℋ′′::IT
    η::IT
    η′::IT
    η′′::IT
    ρ₀ℳ::IT
    # P₀ℳ::IT
end

function Background(par::AbstractCosmoParams{T}; x_grid=-20.0:0.01:0.0, nq=15) where T
    quad_pts, quad_wts =  gausslegendre( nq ) #12 should get 1e-3, 15 conservative
    #Passing the quad pts/wts gets a little busy but eliminates quadgk
    #We may want to fix the quad points to be more/less for bg compared to perts
    #e.g. CLASS uses tolerances of 1e-5 for bg and 1e-3 for perts'

    #println([ρP_0(x2a(x), par,quad_pts,quad_wts) for x in x_grid])
    # println("background T: ", T)

    #FIXME do the tuple juggling to avoid calling quad twice for ρ and P
    ρ₀ℳ_ = spline([ρP_0(x2a(x), par,quad_pts,quad_wts)[1] for x in x_grid], x_grid)
    # P₀ℳ_ = spline([ρP_0(x2a(x), par,quad_pts,quad_wts)[2] for x in x_grid], x_grid)
    ℋ_  = spline([ℋ(x, par,quad_pts,quad_wts) for x in x_grid], x_grid)
    η_   = spline([η(x, par,quad_pts,quad_wts) for x in x_grid], x_grid)
    # println("ℋ_ T: ",typeof(ℋ_))
    return Background(
        T(H₀(par)),
        T(η(0.0, par,quad_pts,quad_wts)),
        T(ρ_crit(par)),
        T(Ω_Λ(par)),

        x_grid,
        # convert(Array{T,1},quad_pts), #explicit call to convert instead of constructor for arrays
        # convert(Array{T,1},quad_wts),
        quad_pts,
        quad_wts,

        ℋ_,
        spline_∂ₓ(ℋ_, x_grid),
        spline_∂ₓ²(ℋ_, x_grid),

        η_,
        spline_∂ₓ(η_, x_grid),
        spline_∂ₓ²(η_, x_grid),
        ρ₀ℳ_,
        # P₀ℳ_,
    )
end
