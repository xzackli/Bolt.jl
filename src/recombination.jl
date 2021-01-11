# Saha Equation
# Useful for high ionization fractions.

# auxillary equations for saha_rhs
n_b(a, par) = par.Ω_b * ρ_crit(par) / (m_H * a^3)
n_H(a, par) = n_b(a, par)  # ignoring helium for now
T_b(a, par) = par.T₀ / a
saha_rhs(a::T, par) where T = (m_e * T_b(a, par) / 2π)^(3/2) / n_H(a, par) *
    exp(-ε₀_H / T_b(a, par))  # rhs of Callin06 eq. 12


function saha_Xₑ(x::T, par::AbstractCosmo{T, DT}) where {T, DT}
    rhs = saha_rhs(x2a(x), par)
    function f!(F, Xe)
        F[1] = Xe[1]^2 / (one(DT) - Xe[1]) - rhs
    end
    res = nlsolve(f!, [0.99 * one(DT)], autodiff = :forward)
    return res.zero[1]
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
Λ_α(a, Xₑ, par) = H(a, par) * (3ε₀_H)^3 / ((8π)^2 * n₁ₛ(a, Xₑ, par))
Cᵣ(a, Xₑ, T_b, par) = (Λ_2s_to_1s + Λ_α(a, Xₑ, par)) / (
    Λ_2s_to_1s + Λ_α(a, Xₑ, par) + β⁽²⁾(T_b))

# RHS of Callin06 eq. 13
function peebles_Xₑ′(Xₑ, par, x)
    a = exp(x)
    T_b_a = BigFloat(T_b(a, par))  # handle overflows by switching to bigfloat
    return float(Cᵣ(a, Xₑ, T_b_a, par) / H(a, par) * (
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
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    return sol  # ode solutions work as interpolator
end


"""
    saha_peebles_recombination(par::AbstractCosmo)

Utility function for generating a decent approximation to Xₑ in ΛCDM recombination,
using the Saha equation until z=1587.4 and then the Peebles equation for the rest.
"""
function saha_peebles_recombination(par::AbstractCosmo{T,DT}) where {T, DT}
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

# visibility functions (OPTIMIZATION OPPORTUNITY) ----------

function τ_integrand_x(x, Xₑ_function, par)
    a = x2a(x)
    return Xₑ_function(x) * n_H(a, par) * σ_T * a / ℋ(a, par)
end
# τ_x(x::Real, Xₑ_function, par) = quadgk(x->τ_integrand_x(x, Xₑ_function, par), x, 0.0)[1]

function τ_x(x::Vector, Xₑ_function, par::AbstractCosmo)
    @assert x[1] > x[2]
    τ_integrands = [-τ_integrand_x(x_, Xₑ_function, par) for x_ in x]
    τ = cumul_integrate(x, τ_integrands)
    return interpolate((reverse(x),),reverse(τ),Gridded(Linear()))
end

function τ̇(x, Xₑ_function, par)
    a = x2a(x)
    return Xₑ_function(x) * n_H(a, par) * a
end

function τ′(x, Xₑ_function, par)
    a = x2a(x)
    return -Xₑ_function(x) * n_H(a, par) * a * σ_T / ℋ(a, par)
end

function g̃(par, Xₑ_function, τ_x_function)
    return x -> -τ′(x, Xₑ_function, par) * exp(-τ_x_function(x))
end

# need τ(x) and g̃(x) for source function
