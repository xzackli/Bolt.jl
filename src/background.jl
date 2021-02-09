
# NOTE: Bolt's background functions are in terms of x ≡ ln(a), the log scale factor

H₀(𝕡::Params) = 𝕡.h * km_s_Mpc_100
ρ_crit(𝕡::Params) = (3 / 8π) * H₀(𝕡)^2 / G_natural  # [eV⁴]
ΩΛ(𝕡::Params) = 1 - (𝕡.Ωr*(1+(7𝕡.Nν/8)*(4/11)^4/3) + 𝕡.Ωb + 𝕡.Ωm)  # dark energy density

#need to account for neutrinos

# Hubble 𝕡ameter ȧ/a in Friedmann background
Hₐ(𝕡::Params, a) = let 
    H₀(𝕡) * √((𝕡.Ωm + 𝕡.Ωb) * a^(-3) + 𝕡.Ωr*(1+(7𝕡.Nν/8)*(4/11)^4/3) * a^(-4) + ΩΛ(𝕡))
end

# conformal time Hubble 𝕡ameter, aH
ℋₐ(𝕡::Params, a) = a * Hₐ(𝕡, a)

# functions in terms of x
H(𝕡::Params, x) = Hₐ(𝕡, x2a(x))
ℋ(𝕡::Params, x) = ℋₐ(𝕡, x2a(x))

# conformal time
η(𝕡::Params, x) = quadgk(a -> 1.0 / (a * ℋₐ(𝕡, a)), 0.0, x2a(x))[1]


# now build a Background with these functions

# a background is 𝕡ametrized on the scalar type T, the interpolator type IT,
# and a type for the grid GT
abstract type AbstractBackground{T, IT<:AbstractInterpolation{T,1}, GT} end

struct Background{T, IT, GT} <: AbstractBackground{T, IT, GT}
    H₀     :: T
    η₀     :: T
    ρ_crit :: T
    ΩΛ     :: T

    x_grid :: GT
    ℋ      :: IT
    ℋ′     :: IT
    ℋ″     :: IT
    η      :: IT
    η′     :: IT
    η″     :: IT
end

function Background(𝕡::Params{T}; x_grid=-20.0:0.01:0.0) where T

    ℋ_ = spline(x_grid, [ℋ(𝕡, x) for x in x_grid])
    η_ = spline(x_grid, [η(𝕡, x) for x in x_grid])

    return Background(
        T(H₀(𝕡)),
        T(η(𝕡, 0.0)),
        T(ρ_crit(𝕡)),
        T(ΩΛ(𝕡)),
        x_grid,

        ℋ_,
        spline_∂ₓ(ℋ_, x_grid),
        spline_∂ₓ²(ℋ_, x_grid),

        η_,
        spline_∂ₓ(η_, x_grid),
        spline_∂ₓ²(η_, x_grid),
    )
end
