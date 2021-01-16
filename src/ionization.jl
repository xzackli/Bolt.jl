
abstract type AbstractIonizationHistory{T, IT<:AbstractInterpolation{T}} end


abstract type IonizationIntegrator end
struct Peebles <: IonizationIntegrator end

struct IonizationHistory{T, IT} <: AbstractIonizationHistory{T, IT}
    Xₑ::IT
    τ::IT
    τ′::IT
    τ′′::IT
    g̃::IT
    g̃′::IT
    g̃′′::IT
end


function IonizationHistory(integrator::Peebles, par::ACP, bg::AB) where
                           {T, ACP<:AbstractCosmoParams{T}, AB<:AbstractBackground}
    x_grid = bg.x_grid
    Xₑ_function = Bolt.saha_peebles_recombination(par)
    τ, τ′ = τ_functions(x_grid, Xₑ_function, par)
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
const T₀ = ustrip(natural(2.725u"K"))  # CMB temperature [K]  # TODO: make this a parameter of the ionization
n_b(a, par) = par.Ω_b * ρ_crit(par) / (m_H * a^3)
n_H(a, par) = n_b(a, par)  # ignoring helium for now
T_b(a, par) = T₀ / a
saha_rhs(a, par) = (m_e * T_b(a, par) / 2π)^(3/2) / n_H(a, par) *
    exp(-ε₀_H / T_b(a, par))  # rhs of Callin06 eq. 12

function saha_Xₑ(x, par::AbstractCosmoParams)
    rhs = saha_rhs(x2a(x), par)
    return  (√(rhs^2 + 4rhs) - rhs) / 2  # solve Xₑ² / (1-Xₑ) = RHS, it's a polynomial
end


"""
    saha_Xₑ(par) where T

Construct an interpolator for mapping scale factor to the ionization fraction
predicted by the Saha equation.

# Arguments:
- `x::T`: log scale factor
- `par`: cosmological parameters structure

# Returns:
- `T`: ionization fraction

# Examples
```julia-repl
julia> f = Bolt.saha_Xₑ(Cosmo()); f(0.00062956)
0.9899574622791693
```
"""
saha_Xₑ(par) = (x -> saha_Xₑ(x, par))


# Peebles Equation
# Use this for Xₑ < 0.99, i.e. z < 1587.4

# recombination parameters for Saha/Peebles
const Λ_2s_to_1s = ustrip(natural(8.227u"s^-1"))  # rate of hydrogen double transition from 2s to 1s
const ε₀_H = ustrip(natural(13.605698u"eV"))  # ionization energy of hydrogen
const m_e = ustrip(natural(float(ElectronMass)))
const m_H = ustrip(natural(float(ProtonMass)))
const α = ustrip(natural(float(FineStructureConstant)))
const σ_T = ustrip(natural(float(ThomsonCrossSection)))

# auxillary equations
ϕ₂(T_b) = 0.448 * log(ε₀_H / T_b)
α⁽²⁾(T_b) = (64π / √(27π)) * (α^2 / m_e^2) * √(ε₀_H / T_b) * ϕ₂(T_b)
β(T_b) = α⁽²⁾(T_b) * (m_e * T_b / (2π))^(3/2) * exp(-ε₀_H / T_b)
β⁽²⁾(T_b) = β(T_b) * exp(3ε₀_H / 4T_b)
n₁ₛ(a, Xₑ, par) = (1 - Xₑ) * n_H(a, par)
Λ_α(a, Xₑ, par) = H_a(a, par) * (3ε₀_H)^3 / ((8π)^2 * n₁ₛ(a, Xₑ, par))
Cᵣ(a, Xₑ, T_b, par) = (Λ_2s_to_1s + Λ_α(a, Xₑ, par)) / (
    Λ_2s_to_1s + Λ_α(a, Xₑ, par) + β⁽²⁾(T_b))

# RHS of Callin06 eq. 13
function peebles_Xₑ′(Xₑ, par, x)
    a = exp(x)
    T_b_a = BigFloat(T_b(a, par))  # handle overflows by switching to bigfloat
    return float(Cᵣ(a, Xₑ, T_b_a, par) / H_a(a, par) * (
        β(T_b_a) * (1 - Xₑ) - n_H(a, par) * α⁽²⁾(T_b_a) * Xₑ^2))
end


"""
    peebles_Xₑ(par, Xₑ₀, x_start, x_end)

Solve the Peebles equation over a span of scale factors, and then
construct an interpolator mapping scale factor to the resulting
ionization fraction.

# Arguments:
- `par`: cosmological parameters
- ` Xₑ₀`: initial ionization fraction
- `x_start`: scale factor to begin integration
- `x_end`: scale factor to end integration

# Returns:
- `generic function`: interpolator for Xₑ(x)
"""
function peebles_Xₑ(par, Xₑ₀, x_start, x_end)
    # set up problem and integrate dXₑ/dx = peebles_Xₑ′
    prob = ODEProblem(peebles_Xₑ′, Xₑ₀, (x_start, x_end), par)
    sol = solve(prob, Tsit5(), reltol=1e-11, abstol=1e-11, dense=true)
    return sol  # ode solutions work as interpolator
end


"""
    saha_peebles_recombination(par::AbstractCosmoParams)

Utility function for generating a decent approximation to Xₑ in ΛCDM recombination,
using the Saha equation until z=1587.4 and then the Peebles equation for the rest.
"""
function saha_peebles_recombination(par::AbstractCosmoParams{T}) where {T}
    z_transition = 1587.4
    x_transition = z2x(z_transition)
    saha_z_grid = 1800:-10:z_transition
    peebles_z_grid = z_transition:-10:100
    early_time_Xₑ = Bolt.saha_Xₑ(par)
    late_time_Xₑ = Bolt.peebles_Xₑ(
        par, early_time_Xₑ(x_transition), x_transition, 0.0)
    Xₑ = x -> (x < x_transition) ? early_time_Xₑ(x) : late_time_Xₑ(x)
    return Xₑ
end

function τ_functions(x, Xₑ_function, par::AbstractCosmoParams)
    @assert x[2] > x[1]  # CONVENTION: x increasing always
    # do a reverse cumulative integrate
    rx = reverse(x)
    τ_primes = [τ′(x_, Xₑ_function, par) for x_ in x]
    τ_integrated = reverse(cumul_integrate(rx, reverse(τ_primes)))

    τ̂ = interpolate((x,),τ_integrated,Gridded(Linear()))
    τ̂′ = interpolate((x,),τ_primes,Gridded(Linear()))
    return τ̂, τ̂′
end

function τ̇(x, Xₑ_function, par)
    a = x2a(x)
    return Xₑ_function(x) * n_H(a, par) * a
end

function τ′(x, Xₑ_function, par)
    a = x2a(x)
    return -Xₑ_function(x) * n_H(a, par) * a * σ_T / ℋ_a(a, par)
end

function g̃_function(τ_x_function, τ′_x_function)
    return x -> -τ′_x_function(x) * exp(-τ_x_function(x))
end
