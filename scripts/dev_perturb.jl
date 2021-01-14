using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using OffsetArrays
using OrdinaryDiffEq
import Bolt: H₀
using ForwardDiff

# utility methods
δ_kron(i, j) = (i == j) ? 1 : 0
∂ₓ(f, x) = ForwardDiff.derivative(f, x)

par = CosmoParams()


k = 340H₀(par)
xgrid = collect(-20.0:0.001:0.0)
zgrid = x2z.(xgrid)
Xₑ = Bolt.saha_peebles_recombination(par)
τ, τ′, τ′′ = Bolt.τ_functions(xgrid, Xₑ, par)
g̃, g̃′ = Bolt.g̃_functions(τ, τ′, τ′′)
g̃′′ = x -> ∂ₓ(g̃′, x)
η = Bolt.η_function(xgrid, par)
ℋ, ℋ′ = x -> Bolt.ℋ(x, par), x -> Bolt.ℋ′(x, par)


TCA_condition(k, ℋₓ, τₓ′) = (abs(k / (ℋₓ * τₓ′)) < 0.1) & (abs(τₓ′) > 10.0)

# NOTE: NO NEUTRINOS 𝒩
function hierarchy!(du, u, p::AbstractCosmoParams, x)
    ℓᵧ, Ω_r, Ω_b, Ω_m = p.ℓᵧ, p.Ω_r, p.Ω_b, p.Ω_m
    H₀² = H₀(p)^2
    ℋₓ, ℋₓ′ = ℋ(x), ℋ′(x)
    τₓ′, τₓ′′ = τ′(x), τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)

    # get array views of photon perturbations
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(view(du, 1:(ℓᵧ+1)), 0:ℓᵧ)
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(view(du, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]

    # metric perturbations
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2])
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b + 4Ω_r * a^(-2) * Θ[0])

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′

    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    if TCA_condition(k, ℋₓ, τₓ′)
        Θ′[2] = 0.0  # could be solved for
        term_3Θ′_vb′ = (
            -((1-R)*τₓ′ + (1+R)*τₓ′′) * (3Θ[1] + v_b) - k * Ψ / ℋₓ
            + (1 - ℋₓ′ / ℋₓ) * (k / ℋₓ) * (-Θ[0] + 2Θ[2]) + k / ℋₓ * (-Θ′[0] + 2Θ′[2])
        ) / ((1+R)*τₓ′ + ℋₓ′ / ℋₓ - 1)
        v_b′ = (-v_b - k * Ψ / ℋₓ + R * (
            term_3Θ′_vb′ + k / ℋₓ * (-Θ[0] + 2Θ[2]) - k / ℋₓ * Ψ
        )) / (1+R)
        Θ′[1] = (term_3Θ′_vb′ - v_b′) / 3
    else
        v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)
        Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    end

    # photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    for ℓ in 2:(ℓᵧ-1)
        Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end
    # polarized photons
    Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓᵧ-1)
        Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end
    # photon hierarchy boundary conditions
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * η(x)) + τₓ′ * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * η(x)) + τₓ′ * Θᵖ[ℓᵧ]

    du[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end


function adiabatic_initial_conditions(par::AbstractCosmoParams{T,DT}, xᵢ) where {T,DT}
    ℓᵧ = par.ℓᵧ
    u = zeros(DT, 2ℓᵧ+7)
    ℋₓ = Bolt.ℋ(xᵢ, par)
    τₓ′ = Bolt.τ′(xᵢ, Xₑ, par)
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ

    # metric and matter perturbations
    Φ = 1.0
    δ = 3Φ / 2
    δ_b = δ
    v = k / (2ℋₓ) * Φ
    v_b = v

    # photon hierarchy
    Θ[0] = Φ / 2
    Θ[1] = -k * Φ / (6ℋₓ)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓᵧ
        Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end

    u[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ, δ, v, δ_b, v_b  # pack back in
    return u
end

xᵢ = log(1e-8)
u₀ = adiabatic_initial_conditions(par, xᵢ)
prob = ODEProblem(hierarchy!, u₀, (xᵢ , 0.0), par)
@time sol = solve(prob, AutoTsit5(Rodas5()), reltol=1e-10)

##
figure()
clf()
sol_x = sol.t
plot(x2a.(sol_x), [abs(sol(x)[3]) for x in sol_x], "-", label=raw"$|\Theta_{\ell=3}|$ photon mode for $k=340H_0$")
xlabel(raw"$a$")
yscale("log")
xscale("log")
# xlim(1e-6, 1e0)
legend()
# ylim(1e-1, 1e8)
gcf()

##
clf()
plot(x2a.(xgrid), [TCA_condition(k, ℋ(x), τ′(x)) for x in xgrid], "-", label="TCA yes no")
axvline(1/1101, ls="dashed", label="recombination")
xscale("log")
xlabel(raw"$a$")
legend()
gcf()

##

##
function source_function(sol, k, x, par)
    u = sol(x)
    u′ = similar(u)
    hierarchy!(u′, u, par, x)

    ℓᵧ = par.ℓᵧ
    H₀² = H₀(par)^2
    ℋₓ, ℋₓ′, ℋₓ′′ = ℋ(x), ℋ′(x), ∂ₓ(ℋ′, x)
    τₓ, τₓ′, τₓ′′ = τ(x), τ′(x), τ′′(x)
    g̃ₓ, g̃ₓ′, g̃ₓ′′ = g̃(x), g̃′(x), g̃′′(x)

    # unpack variables from solution
    Θ = OffsetVector(u[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ = OffsetVector(u[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(u′[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(u′[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]
    Φ′, δ′, v′, δ_b′, v_b′ = u′[(2ℓᵧ+3):(2ℓᵧ+7)]

    a = x2a(x)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * par.Ω_r * Θ[2]
    Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * par.Ω_r * (Θ′[2] - 2 * Θ[2])

    term1 =  g̃ₓ * (Θ[0] + Ψ + Π/4) + exp(-τₓ) * (Ψ′ - Φ′)
    term2 = (-1/k) * (ℋₓ′ * g̃ₓ * v_b + ℋₓ * g̃ₓ′ * v_b + ℋₓ * g̃ₓ * v_b′)
    Π′′ = 2k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * Θ[1] + Θ′[1]) + (3/10) * (τₓ′′ * Π + τₓ′ * Π′) -
        3k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * (Θ[3] + Θᵖ[1] + Θᵖ[3]) + (Θ′[3] + Θᵖ′[1] + Θᵖ′[3]))

    term3 = (3/(4k^2)) * (
        (ℋₓ′^2 + ℋₓ * ℋₓ′′) * g̃ₓ * Π + 3 * ℋₓ * ℋₓ′ * (g̃ₓ′ * Π + g̃ₓ * Π′) +
        ℋₓ^2 * (g̃ₓ′′ * Π + 2g̃ₓ′ * Π′ + g̃ₓ * Π′′)
    )
    return term1 + term2 + term3
end

@time source_function(sol, k, -5.0, par)

##
using SpecialFunctions

S̃(x) = source_function(sol, k, x, par)
η₀ = η(0.0)

clf()
xx = -7.5:0.01:-0.2
ℓ̂ = 100
plot(xx, [S̃(x) * sphericalbesselj(ℓ̂, k*(η₀ - η(x))) for x in xx], "-", lw=0.5)

ylabel(raw"Source function $\times$ bessel")
xlabel(raw"$x$")
gcf()

##
using QuadGK
@time quadgk(x -> S̃(x) * sphericalbesselj(ℓ̂, k*(η₀ - η(x))) , xᵢ, 0.0, rtol=1e-5)

##
