
# NOTE: Bolt's background functions are in terms of x ≡ ln(a), the log scale factor

const ζ = 1.2020569 #Riemann ζ(3) for phase space integrals

H₀(par::AbstractCosmoParams) = par.h * km_s_Mpc_100
ρ_crit(par::AbstractCosmoParams) = (3 / 8π) * H₀(par)^2 / G_natural
function Ω_Λ(par::AbstractCosmoParams)
    #Below can definitely be more streamlined, I am just making it work for now
    return 1 - (par.Ω_r + par.Ω_b + par.Ω_c) 
end

# Hubble parameter ȧ/a in Friedmann background
function H_a(a, par::AbstractCosmoParams)
    return H₀(par) * √(par.Ω_c * a^par.α_c
                        + par.Ω_b * a^(-3)
                        + par.Ω_r* a^(-4)
                        + Ω_Λ(par))
end
# conformal time Hubble parameter, aH
ℋ_a(a, par::AbstractCosmoParams) = a * H_a(a, par)

# functions in terms of x
H(x, par::AbstractCosmoParams) = H_a(x2a(x),par)
ℋ(x, par::AbstractCosmoParams) = ℋ_a(x2a(x), par)

# conformal time
function η(x, par::AbstractCosmoParams,quad_pts,quad_wts)
    logamin,logamax=-13.75,log10(x2a(x))
    Iη(y) = 1.0 / (xq2q(y,logamin,logamax) * ℋ_a(xq2q(y,logamin,logamax), par))/ dxdq(xq2q(y,logamin,logamax),logamin,logamax)
    return sum(Iη.(quad_pts).*quad_wts)
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
    quad_pts::Array{T,1}
    quad_wts::Array{T,1}

    ℋ::IT
    ℋ′::IT
    ℋ′′::IT
    η::IT
    η′::IT
    η′′::IT
end

function Background(par::AbstractCosmoParams{T}; x_grid=-20.0:0.01:0.0, nq=15) where T
    quad_pts, quad_wts =  gausslegendre( nq ) 
    ℋ_  = spline([ℋ(x, par) for x in x_grid], x_grid)
    η_  = spline([η(x, par,quad_pts,quad_wts) for x in x_grid], x_grid)
    return Background(
        T(H₀(par)),
        T(η(0.0, par,quad_pts,quad_wts)),
        T(ρ_crit(par)),
        T(Ω_Λ(par)),

        x_grid,
        convert(Array{T,1},quad_pts), #explicit call to convert instead of constructor for arrays
        convert(Array{T,1},quad_wts),

        ℋ_,
        spline_∂ₓ(ℋ_, x_grid),
        spline_∂ₓ²(ℋ_, x_grid),

        η_,
        spline_∂ₓ(η_, x_grid),
        spline_∂ₓ²(η_, x_grid),
    )
end
