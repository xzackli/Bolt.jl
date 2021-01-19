
abstract type PerturbationIntegrator end
struct BasicNewtonian <: PerturbationIntegrator end

# purely for packing together everything needed to integrate a hierarchy at wavenumber k
struct Hierarchy{T<:Real, Tk<:Real, CP<:AbstractCosmoParams{T},
                 BG<:AbstractBackground, IH<:AbstractIonizationHistory}
    k::Tk
    ℓᵧ::Int  # Boltzmann hierarchy cutoff, i.e. Seljak & Zaldarriaga
    par::CP
    bg::BG
    ih::IH
end

# basic newtonian gauge integrator
function boltsolve(hierarchy::Hierarchy{T}, integrator::BasicNewtonian,
                   ode_alg; reltol=1e-10) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions(hierarchy.par, xᵢ, hierarchy, integrator)
    prob = ODEProblem(basic_newtonian_hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, dense=true)
    return sol
end

# specify a default integrator
boltsolve(hierarchy::Hierarchy, integrator::BasicNewtonian; reltol=1e-11) = boltsolve(
    hierarchy, integrator, Rodas5(); reltol=reltol)


# this integrator comes from Callin+06 and the Dodelson textbook
function basic_newtonian_hierarchy!(du, u, hierarchy::Hierarchy, x)
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_m, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, bg.H₀^2
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)

    # get array views of photon perturbations. this could be sugar'd up probably
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(view(du, 1:(ℓᵧ+1)), 0:ℓᵧ)
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(view(du, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]

    # preliminary variables
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)

    # metric perturbations
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2])
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b + 4Ω_r * a^(-2) * Θ[0])

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)

    # photons and polarized photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    for ℓ in 2:(ℓᵧ-1)
        Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end

    Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓᵧ-1)
        Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end

    # photon boundary conditions: diffusion damping
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θᵖ[ℓᵧ]

    du[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end

# TO IMPROVE: make hierarchy dispatch on solution eltype
function initial_conditions(par::AbstractCosmoParams{T}, xᵢ, hierarchy,
        integrator::BasicNewtonian) where {T}
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    u = zeros(T, 2ℓᵧ+7)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
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

# TODO: this can be extended to any Newtonian gauge integrator if we specify the
# Bardeen potential Ψ and its derivative ψ′
function source_function(du, u, hierarchy, x, par, integrator::BasicNewtonian)
    # compute some quantities
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ′′(x)
    τₓ, τₓ′, τₓ′′ = ih.τ(x), ih.τ′(x), ih.τ′′(x)
    g̃ₓ, g̃ₓ′, g̃ₓ′′ = ih.g̃(x), ih.g̃′(x), ih.g̃′′(x)

    # unpack variables from solution
    Θ = OffsetVector(u[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ = OffsetVector(u[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(du[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(du[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]
    Φ′, δ′, v′, δ_b′, v_b′ = du[(2ℓᵧ+3):(2ℓᵧ+7)]

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
