# Saha Equation
# Useful for high ionization fractions.

# auxillary equations for saha_rhs
n_b(a, par) = par.Ω_b * ρ_crit(par) / (m_H * a^3)
n_H(a, par) = n_b(a, par)  # ignoring helium for now
T_b(a, par) = par.T₀ / a
saha_rhs(a::T, par) where T = (m_e * T_b(a, par) / 2π)^(3/2) / n_H(a, par) *
    exp(-ε₀_H / T_b(a, par))  # rhs of Callin06 eq. 12


function saha_Xₑ(a::T, par::AbstractCosmo{T, DT}) where {T, DT}
    rhs = saha_rhs(a, par)
    function f!(F, x)
        F[1] = x[1]^2 / (one(DT) - x[1]) - rhs
    end
    res = nlsolve(f!, [0.99 * one(DT)], autodiff = :forward)
    return res.zero[1]
end

"""
    saha_Xₑ(par) where T

Compute the ionization fraction at a given scale factor and set of cosmological parameters.

# Arguments:
- `a::T`: scale factor
- `par`: cosmological parameters structure

# Returns:
- `T`: ionization fraction

# Examples
```julia-repl
julia> f = Bolt.saha_Xₑ(Cosmo()); f(0.00062956)
0.9899574622791693
```
"""
saha_Xₑ(par) = (a -> saha_Xₑ(a, par))


# Peebles Equation
# Use this for Xₑ < 0.99, i.e. z < 1587.4

# recombination parameters for Saha/Peebles
const Λ_2s_to_1s = ustrip(natural(8.227u"s^-1"))  # rate of hydrogen double transition from 2s to 1s
const ε₀_H = ustrip(natural(13.605698u"eV"))  # ionization energy of hydrogen
const m_e = ustrip(natural(float(ElectronMass)))
const m_H = ustrip(natural(float(ProtonMass)))
const α = ustrip(natural(float(FineStructureConstant)))

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
    T_b_a = T_b(a, par)
    return Cᵣ(a, Xₑ, T_b_a, par) / H(a, par) * (
        β(T_b_a) * (1 - Xₑ) - n_H(a, par) * α⁽²⁾(T_b_a) * Xₑ^2)
end


"""
    peebles_Xₑ(par, Xₑ₀, a_start, a_end)

Solve the Peebles equation over a span of scale factors, and then
construct an interpolator mapping scale factor to the resulting
ionization fraction.

# Arguments:
- `par`: cosmological parameters
- ` Xₑ₀`: initial ionization fraction
- `a_start`: scale factor to begin integration
- `a_end`: scale factor to end integration

# Returns:
- `generic function`: interpolator for Xₑ(a)
"""
function peebles_Xₑ(par, Xₑ₀, a_start, a_end)
    # integrate dXₑ/dx = peebles_Xₑ′
    prob = ODEProblem(peebles_Xₑ′, Xₑ₀, (log(a_start), log(a_end)), par)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    return a -> sol(log(a))
end


# visibility functions
τ(η, η₀, X⃗ₑ, a⃗) = 1
