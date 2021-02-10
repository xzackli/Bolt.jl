# these types and functions integrate the Boltzmann hierarchy through time

abstract type PerturbationIntegrator end
struct BasicNewtonian <: PerturbationIntegrator end

# a container for everything needed to integrate a hierarchy at wavenumber k
struct Hierarchy{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractParams{T},
                 BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real}
    integrator::PI
    par::CP
    bg::BG
    ih::IH
    k::Tk
    ℓᵧ::Int  # Boltzmann hierarchy cutoff, i.e. Seljak & Zaldarriaga
end
Hierarchy(integrator::PerturbationIntegrator, par::AbstractParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real, ℓᵧ=8) = Hierarchy(integrator, par, bg, ih, k, ℓᵧ)

@⌛ function boltsolve(hierarchy::Hierarchy{T}, ode_alg=Rodas5(); reltol=1e-10) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, saveat=hierarchy.bg.x_grid, dense=false)
    return sol
end

# basic Newtonian gauge: establish the order of perturbative variables in the ODE solve
function unpack(u, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    ℓᵧ = hierarchy.ℓᵧ
    ℓ_ν = 10 #Callin06, for now
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    𝒩 = OffsetVector(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+(ℓ_ν+1))) , 0:ℓ_ν)  # indexed 0 through ℓ_ν
    Φ, δ, v, δ_b, v_b = view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+5)) #getting a little messy...
    return Θ, Θᵖ, 𝒩, Φ, δ, v, δ_b, v_b
end

# BasicNewtonian comes from Callin+06 and the Dodelson textbook (dispatches on hierarchy.integrator)
@⌛ function hierarchy!(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ωr, Ωb, Ωm, Nν, H₀² = par.Ωr, par.Ωb, par.Ωm, par.Nν, bg.H₀^2 #add Nν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ″ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ″(x)
    a = x2a(x)
    R = 4Ωr / (3Ωb * a)
    Ω_ν =  7Nν/8 *(4/11)^(4/3) *Ωr
    ℓ_ν = 10 #again, for now

    Θ, Θᵖ, 𝒩, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′,_, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6

    # metric perturbations
    #Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ωr * Θ[2])
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ωr * Θ[2] + Ω_ν * 𝒩[2]) #add rel quadrupole
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ωm * a^(-1) * δ + Ωb * a^(-1) * δ_b + 4Ωr * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0]) #add rel monopole on this line

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)

    # relativistic neutrinos (massless)
    𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
    𝒩′[1] = k/(3ℋₓ)*𝒩[0] - 2*k/(3ℋₓ)*𝒩[2] + k/(3ℋₓ)*Ψ
    for ℓ in 2:(ℓ_ν-1) #ℓ_ν same as ℓᵧ for massless nu for now
        𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) *( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1])
    end
    #truncation
    𝒩′[ℓ_ν] =  k / ℋₓ  * 𝒩[ℓ_ν-1] - (ℓ_ν+1)/(ℋₓ *ηₓ) *𝒩[ℓ_ν]#Callin 06

    # WIP: massive nu
    # # neutrinos (massive, MB 57) - change convention
    # #units not yet right
    # Ψ_ν′[0] = -q*k / ϵ * Ψ_ν[1] - Φ′ *dnlnf0dlnq #FIXME dnln, def Psi,IC, Einstein int
    # Ψ_ν′[1] = q*k/(3 ϵ)*(𝒩[0] - 2*Ψ_ν[2]) - This k*ϵ / (3*q)*Ψ *dnlnf0dlnq
    # I think can't mutate the u variable...
    # Ψ_ν[ℓ_νm] = (2*ℓ_νm+1)*ϵ/(q*k*ηₓ)*Ψ_ν[ℓ_νm] - Ψ_ν[ℓ_νm-1] #truncation of MB (51)
    # for ℓ in 2:(ℓ_νm-1) #ℓ_νm should be smaller than massless case
    #     Ψ_ν′[ℓ] =  q*k / ((2ℓ+1) * ϵ) *( ℓ*Ψ_ν[ℓ-1] - (ℓ+1)*Ψ_ν[ℓ+1])
    # end

    # photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
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

    # photon boundary conditions: diffusion damping
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θᵖ[ℓᵧ]

    du[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end

# BasicNewtonian Integrator (dispatches on hierarchy.integrator)
@⌛ function initial_conditions(xᵢ, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    ℓ_ν = 10 #again, for now
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ″ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ″(xᵢ)
    Θ, Θᵖ, 𝒩, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2


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

    # neutrino hierarchy
    # for now we assume xᵢ is before neutrinos decouple
    f_ν = 1/(1 + 1/(7par.Nν/8 *(4/11)^(4/3)))
    𝒩[0] = Θ[0]
    𝒩[1] = Θ[1]
    𝒩[2] = - (k^2 *aᵢ²*Φ) / (12H₀²) * 1 / (1 + 5f_ν/2) #Callin06 (71)
    for ℓ in 3:ℓ_ν
        𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #approximation of Callin06 (72)
    end

    #WIP: massive nu
    #FIXME: nonrelativistic transition for massive species, needs to go in bg
    #^this will have to wait for m_ν to be added to pars
    #below notation is not right yet, starting from MB
    #x_nr = m_ν/5.3e-4 -1 #m_ν in eV (PDG26-pg3)
    # same as photons for 0,1
    # # massive #FIXME ingegrate the q moments, get dlnf0dlnq,define ϵ
    # σ_ν= (k*ηₓ)^2 *Ψ/ 15
    # Ψ_ν[0] = -δ_ν *dlnf0dlnq
    # Ψ_ν[1] = -ϵ/(3*q*k) θ_ν *dlnf0dlnq  #change θ to -k/ℋₓ v
    # Ψ_ν[2] = -σ_ν/2  *dlnf0dlnq
    #ignore ℓ>2, small

    u[(2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end

# TODO: this could be extended to any Newtonian gauge integrator if we specify the
# Bardeen potential Ψ and its derivative ψ′ for an integrator, or we saved them
function source_function(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute some quantities
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ″ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ″(x)
    τₓ, τₓ′, τₓ″ = ih.τ(x), ih.τ′(x), ih.τ″(x)
    g̃ₓ, g̃ₓ′, g̃ₓ″ = ih.g̃(x), ih.g̃′(x), ih.g̃″(x)
    a = x2a(x)

    Θ, Θᵖ, 𝒩, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ′, Θᵖ′, 𝒩′, Φ′, δ′, v′, δ_b′, v_b′ = unpack(du, hierarchy)

    # recalulate these since we didn't save them
    Ψ = -Φ - 12H₀² / k^2 / a^2 * par.Ωr * Θ[2]
    Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * par.Ωr * (Θ′[2] - 2 * Θ[2])
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]

    term1 =  g̃ₓ * (Θ[0] + Ψ + Π/4) + exp(-τₓ) * (Ψ′ - Φ′)
    term2 = (-1/k) * (ℋₓ′ * g̃ₓ * v_b + ℋₓ * g̃ₓ′ * v_b + ℋₓ * g̃ₓ * v_b′)
    Π″ = 2k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * Θ[1] + Θ′[1]) + (3/10) * (τₓ″ * Π + τₓ′ * Π′) -
        3k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * (Θ[3] + Θᵖ[1] + Θᵖ[3]) + (Θ′[3] + Θᵖ′[1] + Θᵖ′[3]))
    term3 = (3/(4k^2)) * (
        (ℋₓ′^2 + ℋₓ * ℋₓ″) * g̃ₓ * Π + 3 * ℋₓ * ℋₓ′ * (g̃ₓ′ * Π + g̃ₓ * Π′) +
        ℋₓ^2 * (g̃ₓ″ * Π + 2g̃ₓ′ * Π′ + g̃ₓ * Π″))
    return term1 + term2 + term3
end
