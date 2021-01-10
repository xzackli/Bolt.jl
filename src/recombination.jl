
# Saha Equation
# Useful for high ionization fractions.

# auxillary equations for Saha
n_b(a, par) = par.Ω_b * par.ρ_c / (par.m_H * a^3)
T_b(a, par) = par.T₀ / a
saha_rhs(a::T, par) where T = (par.m_e * T_b(a, par) / 2π)^(3/2) / n_b(a, par) *
    exp(-par.ε₀_H / T_b(a, par))

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
1-element Array{Float64,1}:
 0.9899574622791693
```
"""
function saha_Xₑ(a::T, par) where {T <: Real}
    rhs = saha_rhs(a, par)
    function f!(F, x)
        F[1] = x[1]^2 / (one(T) - x[1]) - rhs
    end
    res = nlsolve(f!, [0.99], autodiff = :forward)
    return res.zero[1]
end
saha_Xₑ(a⃗::Vector{T}, par) where T = T[saha_Xₑ(a, par) for a in a⃗]

# Peebles Equation
# Use this for Xₑ < 0.99, i.e. z < 1587.4

# auxillary equations
ϕ₂(T_b, par) = 0.448 * log(par.ε₀_H / T_b)
α⁽²⁾(T_b, par) = (64π / √(27π)) * (par.α^2 / par.m_e^2) * √(par.ε₀_H / T_b) * ϕ₂(T_b, par)
β(T_b, par) = α⁽²⁾(T_b, par) * (par.mₑ * T_b / (2π))^(3/2) * exp(-par.ε₀_H / T_b)
β⁽²⁾(T_b, par) = β(T_b, par) * exp(3par.ε₀_H / 4T_b)
# n₁ₛ(Xₑ::T) where T = (one(T) - Xₑ) * n_H

function peebles_Xₑ()
end
