# recombination parameters for Saha/Peebles
const Λ_2s_to_1s = ustrip(natural(8.227u"s^-1"))  # rate of hydrogen double transition from 2s to 1s
const ε₀_H = ustrip(natural(13.605698u"eV"))  # ionization energy of hydrogen
const m_e = ustrip(natural(float(ElectronMass)))
const m_H = ustrip(natural(float(ProtonMass)))
const α = ustrip(natural(float(FineStructureConstant)))

# Saha Equation
# Useful for high ionization fractions.

# auxillary equations for saha_rhs
n_b(a, par) = par.Ω_b * par.ρ_c / (m_H * a^3)
n_H(a, par) = n_b(a, par)  # ignoring helium for now
T_b(a, par) = par.T₀ / a
saha_rhs(a::T, par) where T = (m_e * T_b(a, par) / 2π)^(3/2) / n_H(a, par) *
    exp(-ε₀_H / T_b(a, par))  # rhs of Callin06 eq. 12

"""
    saha_Xₑ(a::T, par) where T

Compute the ionization fraction at a given scale factor and set of cosmological parameters.

# Arguments:
- `a::T`: scale factor
- `par`: cosmological parameters structure

# Returns:
- `T`: ionization fraction

# Examples
```julia-repl
julia> Bolt.saha_Xₑ(0.00062956, Cosmo())
0.9899574622791693
```
"""
function saha_Xₑ(a::T, par) where {T <: Real}
    rhs = saha_rhs(a, par)
    function f!(F, x)
        F[1] = x[1]^2 / (one(T) - x[1]) - rhs
    end
    res = nlsolve(f!, [0.99 * one(T)], autodiff = :forward)
    return res.zero[1]
end
saha_Xₑ(a⃗::Vector{T}, par) where T = T[saha_Xₑ(a, par) for a in a⃗]

# Peebles Equation
# Use this for Xₑ < 0.99, i.e. z < 1587.4

# auxillary equations
ϕ₂(T_b) = 0.448 * log(ε₀_H / T_b)
α⁽²⁾(T_b) = (64π / √(27π)) * (α^2 / m_e^2) * √(ε₀_H / T_b) * ϕ₂(T_b)
β(T_b) = α⁽²⁾(T_b) * (m_e * T_b / (2π))^(3/2) * exp(-ε₀_H / T_b)
β⁽²⁾(T_b) = β(T_b) * exp(ε₀_H / 4T_b)
n₁ₛ(a, Xₑ, par) = (1 - Xₑ) * n_H(a, par)
Λ_α(a, Xₑ, par) = H(a, par) * (3ε₀_H)^3 / ((8π)^2 * n₁ₛ(a, Xₑ, par))
Cᵣ(a, Xₑ, T_b, par) = (Λ_2s_to_1s + Λ_α(a, Xₑ, par)) / (
    Λ_2s_to_1s + Λ_α(a, Xₑ, par) + β⁽²⁾(T_b))

# RHS of Callin06 eq. 13
function peebles_rhs(Xₑ, a, par)
    T_b_a = T_b(a, par)
    return Cᵣ(a, Xₑ, T_b_a, par) / H(a, par) * (
        β(T_b_a) * (1 - Xₑ) - n_H(a, par) * α⁽²⁾(T_b_a) * Xₑ^2)
end

function peebles_Xₑ(par, Xₑ₀, a_start, a_end) where T
    function peebles_eq(u,p,t)
        # return Bolt.peebles_rhs(u, t, par)
        return 0.1
    end
    prob = ODEProblem(peebles_eq, Xₑ₀, (a_start, a_end))
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    return sol
end
