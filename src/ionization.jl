
abstract type AbstractIonizationHistory{T, IT<:AbstractInterpolation{T}} end

abstract type IonizationIntegrator end
struct Peebles <: IonizationIntegrator end

struct IonizationHistory{T, IT} <: AbstractIonizationHistory{T, IT}
    Xₑ :: IT
    τ  :: IT
    τ′ :: IT
    τ″ :: IT
    g̃  :: IT
    g̃′ :: IT
    g̃″ :: IT
end

# this Peebles history comes from Callin+06, peep the plots from examples/
# which match that paper perfectly
function IonizationHistory(
    integrator :: Peebles, 
    𝕡 :: Params{T}, 
    bg :: AbstractBackground
) where {T}

    x_grid = bg.x_grid
    Xₑ_function = Bolt.saha_peebles_recombination(𝕡)
    τ, τ′ = τ_functions(𝕡, x_grid, Xₑ_function)
    g̃ = g̃_function(τ, τ′)

    Xₑ_ = spline(x_grid, Xₑ_function.(x_grid))
    τ_ = spline(x_grid, τ.(x_grid))
    g̃_ = spline(x_grid, g̃.(x_grid))
    IT = typeof(Xₑ_)

    # TO FIX, WHY DOES THIS CONSTRUCTOR REQUIRE {I, IT}???
    return IonizationHistory{T, IT}(
        Xₑ_,
        τ_,
        spline_∂ₓ(τ_, x_grid),
        spline_∂ₓ²(τ_, x_grid),
        g̃_,
        spline_∂ₓ(g̃_, x_grid),
        spline_∂ₓ²(g̃_, x_grid),
    )
    
end

# Saha Equation
# Useful for high ionization fractions.

# auxillary equations for saha_rhs
n_b(𝕡, a) = 𝕡.Ωb * ρ_crit(𝕡) / (m_H * a^3)
n_H(𝕡, a) = n_b(𝕡, a)  # ignoring helium for now
Tb(𝕡, a) = T₀ / a
saha_rhs(𝕡, a) = (m_e * Tb(𝕡, a) / 2π)^(3/2) / n_H(𝕡, a) *
    exp(-ε₀_H / Tb(𝕡, a))  # rhs of Callin06 eq. 12

function saha_Xₑ(𝕡::Params, x)
    rhs = saha_rhs(𝕡, x2a(x))
    return  (√(rhs^2 + 4rhs) - rhs) / 2  # solve Xₑ² / (1-Xₑ) = RHS, it's a polynomial
end
saha_Xₑ(𝕡) = (x -> saha_Xₑ(𝕡, x))


# Peebles Equation
# Use this for Xₑ < 0.99, i.e. z < 1587.4


# auxillary equations
ϕ₂(Tb) = 0.448 * log(ε₀_H / Tb)
α⁽²⁾(Tb) = (64π / √(27π)) * (α^2 / m_e^2) * √(ε₀_H / Tb) * ϕ₂(Tb)
β(Tb) = α⁽²⁾(Tb) * (m_e * Tb / (2π))^(3/2) * exp(-ε₀_H / Tb)
β⁽²⁾(Tb) = β(Tb) * exp(3ε₀_H / 4Tb)
n₁ₛ(𝕡, a, Xₑ) = (1 - Xₑ) * n_H(𝕡, a)
Λ_α(𝕡, a, Xₑ) = Hₐ(𝕡, a) * (3ε₀_H)^3 / ((8π)^2 * n₁ₛ(𝕡, a, Xₑ))
Cᵣ(𝕡, a, Xₑ, Tb) = (Λ_2s_to_1s + Λ_α(𝕡, a, Xₑ)) / (
    Λ_2s_to_1s + Λ_α(𝕡, a, Xₑ) + β⁽²⁾(Tb))

# RHS of Callin06 eq. 13
function peebles_Xₑ′(Xₑ, 𝕡, x)
    a = exp(x)
    Tb_a = BigFloat(Tb(𝕡, a))  # handle overflows by switching to bigfloat
    return float(Cᵣ(𝕡, a, Xₑ, Tb_a) / Hₐ(𝕡, a) * (
        β(Tb_a) * (1 - Xₑ) - n_H(𝕡, a) * α⁽²⁾(Tb_a) * Xₑ^2))
end


"""
    peebles_Xₑ(𝕡, Xₑ₀, x_start, x_end)

Solve the Peebles equation over a span of scale factors, and then
construct an interpolator mapping scale factor to the resulting
ionization fraction.

# Arguments:
- `𝕡`: cosmological 𝕡ameters
- ` Xₑ₀`: initial ionization fraction
- `x_start`: scale factor to begin integration
- `x_end`: scale factor to end integration

# Returns:
- `generic function`: interpolator for Xₑ(x)
"""
function peebles_Xₑ(𝕡, Xₑ₀, x_start, x_end)
    # set up problem and integrate dXₑ/dx = peebles_Xₑ′
    prob = ODEProblem(peebles_Xₑ′, Xₑ₀, (x_start, x_end), 𝕡)
    sol = solve(prob, Tsit5(), reltol=1e-11, abstol=1e-11, dense=true)
    return sol  # ode solutions work as interpolator
end


"""
    saha_peebles_recombination(𝕡::AbstractCosmoParams)

Utility function for generating a decent approximation to Xₑ in ΛCDM recombination,
using the Saha equation until z=1587.4 and then the Peebles equation for the rest.
"""
function saha_peebles_recombination(𝕡::Params{T}) where {T}
    z_transition = 1587.4
    x_transition = z2x(z_transition)
    saha_z_grid = 1800:-10:z_transition
    peebles_z_grid = z_transition:-10:100
    early_time_Xₑ = Bolt.saha_Xₑ(𝕡)
    late_time_Xₑ = Bolt.peebles_Xₑ(
        𝕡, early_time_Xₑ(x_transition), x_transition, 0.0)
    Xₑ = x -> (x < x_transition) ? early_time_Xₑ(x) : late_time_Xₑ(x)
    return Xₑ
end

function τ_functions(𝕡, x, Xₑ_function)
    @assert x[2] > x[1]  # CONVENTION: x increasing always
    # do a reverse cumulative integrate
    rx = reverse(x)
    τ_primes = [τ′(𝕡, x_, Xₑ_function) for x_ in x]
    τ_integrated = reverse(cumul_integrate(rx, reverse(τ_primes)))

    τ̂ = interpolate((x,),τ_integrated,Gridded(Linear()))
    τ̂′ = interpolate((x,),τ_primes,Gridded(Linear()))
    return τ̂, τ̂′
end

function τ̇(𝕡, x, Xₑ_function)
    a = x2a(x)
    return Xₑ_function(x) * n_H(𝕡, a) * a
end

function τ′(𝕡, x, Xₑ_function)
    a = x2a(x)
    return -Xₑ_function(x) * n_H(𝕡, a) * a * σ_T / ℋₐ(𝕡, a)
end

function g̃_function(τ_x_function, τ′_x_function)
    return x -> -τ′_x_function(x) * exp(-τ_x_function(x))
end
