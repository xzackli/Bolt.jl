
# derived quantities (I've chosen natural units, possibly the wrong choice)
const km_s_Mpc_100 = ustrip(natural(100.0u"km/s/Mpc"))  # [eV]
const G_natural = ustrip(natural(float(NewtonianConstantOfGravitation)))

H₀(par::AbstractCosmoParams) = par.h * km_s_Mpc_100
ρ_crit(par::AbstractCosmoParams) = (3 / 8π) * H₀(par)^2 / G_natural  # [eV⁴]
Ω_Λ(par::AbstractCosmoParams) = 1 - (par.Ω_r + par.Ω_b + par.Ω_m)  # dark energy density

# Hubble parameter ȧ/a in Friedmann background
H_a(a, par::AbstractCosmoParams) = H₀(par) * √((par.Ω_m + par.Ω_b) * a^(-3) + par.Ω_r * a^(-4) + Ω_Λ(par))
# conformal time Hubble parameter, aH
ℋ_a(a, par::AbstractCosmoParams) = a * H_a(a, par)

# functions in terms of x
H(x, par::AbstractCosmoParams) = H_a(x2a(x),par)
ℋ(x, par::AbstractCosmoParams) = ℋ_a(x2a(x), par)

## hand-written ℋ derivative for testing
# function ℋ′(x, par::AbstractCosmoParams)
#     a = x2a(x)
#     return -H₀(par) * (2par.Ω_r + (par.Ω_b + par.Ω_m) * a - 2Ω_Λ(par) * a^4) /
#         (2 * a * √(par.Ω_r + (par.Ω_b + par.Ω_m) * a + Ω_Λ(par) * a^4))
# end

# conformal time
function η(x, par::AbstractCosmoParams)
    return quadgk(a -> 1.0 / (a * ℋ_a(a, par)), 0.0, x2a(x))[1]
end


# now build a Background with these functions

# a background is parametrized on the scalar type T, the interpolator type IT,
# and a type for the grid GT
abstract type AbstractBackground{T, IT<:AbstractInterpolation{T,1}, GT} end

struct Background{T, IT, GT} <: AbstractBackground{T, IT, GT}
    H₀::T
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
