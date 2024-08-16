## Saha Equation
## Useful for high ionization fractions.

# auxillary equations for saha_rhs
const PeeblesT₀ = ustrip(natural(2.725u"K"))  # CMB temperature [K]  # TODO: make this a parameter of the ionization
n_b(a, par) = par.Ω_b * ρ_crit(par) / (m_H * a^3)
n_H(a, par) = n_b(a, par) *(1-par.Y_p) #Adding Helium!
saha_T_b(a, par) = PeeblesT₀ / a #j why does this take par?
saha_rhs(a, par) = (m_e * saha_T_b(a, par) / 2π)^(3/2) / n_H(a, par) *
    exp(-ε₀_H / saha_T_b(a, par))  # rhs of Callin06 eq. 12

function saha_Xₑ(x, par::AbstractCosmoParams)
    rhs = saha_rhs(x2a(x), par)
    return  (√(rhs^2 + 4rhs) - rhs) / 2  # solve Xₑ² / (1-Xₑ) = RHS, it's a polynomial
end
saha_Xₑ(par) = (x -> saha_Xₑ(x, par))


## Peebles Equation
## Use this for Xₑ < 0.99, i.e. z < 1587.4

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
#Problem is here \/ since Lyα rate is given by redshifting out of line need H
Λ_α(a, Xₑ, par) = oldH_a(a, par) * (3ε₀_H)^3 / ((8π)^2 * n₁ₛ(a, Xₑ, par))
new_Λ_α(a, Xₑ, par, ℋ_function) = ℋ_function(log(a)) * (3ε₀_H)^3 / ((8π)^2 * n₁ₛ(a, Xₑ, par))
Cᵣ(a, Xₑ, T_b, par) = (Λ_2s_to_1s + Λ_α(a, Xₑ, par)) / (
    Λ_2s_to_1s + Λ_α(a, Xₑ, par) + β⁽²⁾(T_b))
new_Cᵣ(a, Xₑ, T_b, par,ℋ_function) = (Λ_2s_to_1s + new_Λ_α(a, Xₑ, par,ℋ_function)) / (
    Λ_2s_to_1s + new_Λ_α(a, Xₑ, par,ℋ_function) + β⁽²⁾(T_b))

# RHS of Callin06 eq. 13
function peebles_Xₑ′(Xₑ, par, x)
    a = exp(x)
    T_b_a = BigFloat(saha_T_b(a, par))  # handle overflows by switching to bigfloat
    return float(Cᵣ(a, Xₑ, T_b_a, par) / oldH_a(a, par) * (
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

function oldτ_functions(x, Xₑ_function, par::AbstractCosmoParams)
    @assert x[2] > x[1]  # CONVENTION: x increasing always
    # do a reverse cumulative integrate
    rx = reverse(x)
    τ_primes = [oldτ′(x_, Xₑ_function, par) for x_ in x]
    τ_integrated = reverse(cumul_integrate(rx, reverse(τ_primes)))

    τ̂ = interpolate((x,),τ_integrated,Gridded(Linear()))
    τ̂′ = interpolate((x,),τ_primes,Gridded(Linear()))
    return τ̂, τ̂′
end

function τ_functions(x, Xₑ_function, par::AbstractCosmoParams,ℋ_function)
    @assert x[2] > x[1]  # CONVENTION: x increasing always
    # do a reverse cumulative integrate
    rx = reverse(x)
    τ_primes = [τ′(x_, Xₑ_function, par, ℋ_function) for x_ in x]
    τ_integrated = reverse(cumul_integrate(rx, reverse(τ_primes)))

    τ̂ = interpolate((x,),τ_integrated,Gridded(Linear()))
    τ̂′ = interpolate((x,),τ_primes,Gridded(Linear()))
    return τ̂, τ̂′
end

function τ̇(x, Xₑ_function, par)
    a = x2a(x)
    return Xₑ_function(x) * n_H(a, par) * a
end

function τ′(x, Xₑ_function, par, ℋ_function)
    a = x2a(x)
    return -Xₑ_function(x) * n_H(a, par) * a * σ_T / ℋ_function(x)
end
function oldτ′(x, Xₑ_function, par)
    a = x2a(x)
    return -Xₑ_function(x) * n_H(a, par) * a * σ_T / (a*oldH_a(a,par))
end

function g̃_function(τ_x_function, τ′_x_function)
    return x -> -τ′_x_function(x) * exp(-τ_x_function(x))
end


"""Convenience function to create an ionisation history from some tables"""
function customion(par, bg, Xₑ_function, Tmat_function, csb²_function)
     x_grid = bg.x_grid
     τ, τ′ = Bolt.τ_functions(x_grid, Xₑ_function, par, bg.ℋ)
     g̃ = Bolt.g̃_function(τ, τ′)
     spline, spline_∂ₓ, spline_∂ₓ² = Bolt.spline, Bolt.spline_∂ₓ, Bolt.spline_∂ₓ²
     Xₑ_ = spline(Xₑ_function.(x_grid), x_grid)
     τ_ = spline(τ.(x_grid), x_grid)
     g̃_ = spline(g̃.(x_grid), x_grid)
     Tmat_ = spline(Tmat_function.(x_grid), x_grid)
     csb²_ = spline(csb²_function.(x_grid), x_grid)
 
     return IonizationHistory(
           (τ(0.)),
         Xₑ_,
         τ_,
         spline_∂ₓ(τ_, x_grid),
         spline_∂ₓ²(τ_, x_grid),
         g̃_,
         spline_∂ₓ(g̃_, x_grid),
         spline_∂ₓ²(g̃_, x_grid),
         Tmat_,
           csb²_,
    )
end
