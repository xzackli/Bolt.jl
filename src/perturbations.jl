# these types and functions integrate the Boltzmann hierarchy through time

abstract type PerturbationIntegrator end
struct BasicNewtonian <: PerturbationIntegrator end

# a container for everything needed to integrate a hierarchy at wavenumber k
struct Hierarchy{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
                 BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real}
    integrator::PI
    par::CP
    bg::BG
    ih::IH
    k::Tk
    ℓᵧ::Int  # Boltzmann hierarchy cutoff, i.e. Seljak & Zaldarriaga
end

Hierarchy(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real, ℓᵧ=8) = Hierarchy(integrator, par, bg, ih, k, ℓᵧ)


struct Hierarchy_nn{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
        BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real,
        # S<:Real,AT<:AbstractArray{S,1}
        }
integrator::PI
par::CP
bg::BG
ih::IH
k::Tk
# p::AT #the 
p::AbstractArray #the 
ℓᵧ::Int  # Boltzmann hierarchy cutoff, i.e. Seljak & Zaldarriaga
end

Hierarchy_nn(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
ih::AbstractIonizationHistory, k::Real, p::AbstractArray,ℓᵧ=8
) = Hierarchy_nn(integrator, par, bg, ih, k, p,ℓᵧ)


struct Hierarchy_spl{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
    BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real,
    S<:Real,IT<:AbstractInterpolation{S,1}}
integrator::PI
par::CP
bg::BG
ih::IH
k::Tk
spl1::IT  
# spl2::IT 
ℓᵧ::Int  # Boltzmann hierarchy cutoff, i.e. Seljak & Zaldarriaga
end

Hierarchy_spl(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
ih::AbstractIonizationHistory, k::Real, 
spl1::AbstractInterpolation,spl2::AbstractInterpolation,ℓᵧ=8
) = Hierarchy_spl(integrator, par, bg, ih, k, spl1,spl2,ℓᵧ)


function boltsolve(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6, abstol=1e-6) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, abstol=abstol,
                # saveat=hierarchy.bg.x_grid, dense=false,  # don't save a grid
                )
    return sol
end

function boltsolve_nn(hierarchy::Hierarchy_nn{T}, ode_alg=KenCarp4(); reltol=1e-6, abstol=1e-6) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions_nn(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy_nn!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, abstol=abstol,
                saveat=hierarchy.bg.x_grid,  
                )
    return sol
end

function boltsolve_spl(hierarchy::Hierarchy_spl{T}, saveat::Array{T,1},ode_alg=KenCarp4(); reltol=1e-6, abstol=1e-6) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions_spl(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy_spl!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, abstol=abstol,
                saveat=saveat,  
                )
    return sol
end

# basic Newtonian gauge: establish the order of perturbative variables in the ODE solve
function unpack(u, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    ℓᵧ = hierarchy.ℓᵧ
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Φ, δ, v, δ_b, v_b = view(u, (2(ℓᵧ+1)+1):(2(ℓᵧ+1)+5)) #getting a little messy...
    return Θ, Θᵖ, Φ, δ, v, δ_b, v_b
end

function unpack_nn(u, hierarchy::Hierarchy_nn{T, BasicNewtonian}) where T
    ℓᵧ = hierarchy.ℓᵧ
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Φ, δ, v, δ_b, v_b = view(u, (2(ℓᵧ+1)+1):(2(ℓᵧ+1)+5)) #getting a little messy...
    return Θ, Θᵖ, Φ, δ,  v, δ_b, v_b
end

function unpack_spl(u, hierarchy::Hierarchy_spl{T, BasicNewtonian}) where T
    ℓᵧ = hierarchy.ℓᵧ
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Φ, δ, v, δ_b, v_b = view(u, (2(ℓᵧ+1)+1):(2(ℓᵧ+1)+5)) #getting a little messy...
    return Θ, Θᵖ, Φ, δ,  v, δ_b, v_b
end

# BasicNewtonian comes from Callin+06 and the Dodelson textbook (dispatches on hierarchy.integrator)
function hierarchy!(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)

    α_c = par.α_c

    Θ, Θᵖ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, _, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6


    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ 
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        )

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)
    

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
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[ℓᵧ]
    #END RSA

    du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end

function get_nn(m,d_in,d_out;rng=Xoshiro(123))
    NN = SimpleChain(static(d_in),
        TurboDense{true}(tanh, m), 
        TurboDense{true}(tanh, m), 
        TurboDense{false}(identity, d_out) #have not tested non-scalar output
    );
    p = SimpleChains.init_params(NN;rng); 
    G = SimpleChains.alloc_threaded_grad(NN);
    return NN
end

function hierarchy_nn!(du, u, hierarchy::Hierarchy_nn{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)
    pertlen=2(ℓᵧ+1)+5

    α_c = par.α_c

    # Θ, Θᵖ, Φ, δ_c, σ_c,δ_b, v_b = unpack_nn(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ, Θᵖ, Φ, δ_c, v_c,δ_b, v_b = unpack_nn(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)

    Θ′, Θᵖ′, _, _, _, _, _ = unpack_nn(du, hierarchy)  # will be sweetened by .. syntax in 1.6


    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                #   + Ω_c * a^(4+α_c) * σ_c 
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ_c
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        )

    # matter
    # THIS DOESN'T WORK!!!
    NN = get_nn(4,pertlen+2,32)
    nnin = hcat([u...,k,x])
    u′ = NN(nnin,hierarchy.p)
    δ′ = u′[1] #k / ℋₓ * v - 3Φ′
    # v′ = -v - k / ℋₓ * Ψ
    v′ = u′[2]
    # σ′ = NN₂

    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)
    

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
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[ℓᵧ]
    #END RSA

    du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    # du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, σ′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end


function hierarchy_spl!(du, u, hierarchy::Hierarchy_spl{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)
    pertlen=2(ℓᵧ+1)+5

    α_c = par.α_c

    # Θ, Θᵖ, Φ, δ_c, σ_c,δ_b, v_b = unpack_nn(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ, Θᵖ, Φ, δ_c, v_c,δ_b, v_b = unpack_spl(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)

    Θ′, Θᵖ′, _, _, _, _, _ = unpack_spl(du, hierarchy)  # will be sweetened by .. syntax in 1.6


    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                #   + Ω_c * a^(4+α_c) * σ_c 
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ_c
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        )

    # matter
    δ′ = hierarchy.spl1(x)
    v′ = -v - k / ℋₓ * Ψ
    # v′ = hierarchy.spl2(x)
    # σ′ = NN₂
    #FIXME we don't actually need to learn v in the ODE, just a scalar at z=0...
    #Well this depends on the assumption we make I think...

    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)
    

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
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[ℓᵧ]
    #END RSA

    du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    # du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, σ′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end


# BasicNewtonian Integrator (dispatches on hierarchy.integrator)

function initial_conditions(xᵢ, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    u = zeros(T, 2(ℓᵧ+1)+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    #These get a 3/3 since massive neutrinos behave as massless at time of ICs
    f_ν = 1/(1 + 1/(7*(3/3)*3.046/8 *(4/11)^(4/3))) # we need to actually keep this

    # metric and matter perturbations
    ℛ = 1.0  # set curvature perturbation to 1
    Φ = (4f_ν + 10) / (4f_ν + 15) * ℛ  # for a mode outside the horizon in radiation era
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) ) * Φ

    #trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓᵧ
        Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end

    δ = 3/4 *(4Θ[0]) #the 4 converts δγ_MB -> Dodelson convention
    δ_b = δ
    #we have that Θc = Θb = Θγ = Θν, but need to convert Θ = - k v (i absorbed in v)
    v = -3k*Θ[1]
    v_b = v

    u[2(ℓᵧ+1)+1:(2(ℓᵧ+1)+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end


function initial_conditions_nn(xᵢ, hierarchy::Hierarchy_nn{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    u = zeros(T, 2(ℓᵧ+1)+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    # Θ, Θᵖ, Φ, δ, σ, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ, Θᵖ, Φ, δ, v, δ_b, v_b = unpack_nn(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    #These get a 3/3 since massive neutrinos behave as massless at time of ICs
    f_ν = 1/(1 + 1/(7*(3/3)*3.046/8 *(4/11)^(4/3))) # we need to actually keep this
    α_c = par.α_c

    # metric and matter perturbations
    ℛ = 1.0  # set curvature perturbation to 1
    Φ = (4f_ν + 10) / (4f_ν + 15) * ℛ  # for a mode outside the horizon in radiation era
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) ) * Φ

    #trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓᵧ
        Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end


    # δ = 3/4 *(4Θ[0]) #the 4 converts δγ_MB -> Dodelson convention
    δ = -α_c*Θ[0] # this is general enough to allow this to be any species
    δ_b = δ
    #we have that Θc = Θb = Θγ = Θν, but need to convert Θ = - k v (i absorbed in v)
    # v = -3k*Θ[1]
    v = α_c*k*Θ[1]
    v_b = -3k*Θ[1]
    # σ = 0.0 # this is an actual physical assumption - that DM has no anisotropic stress in Read
    #^This is strong but oh well, we can revisit making it a free parameter later

    # u[2(ℓᵧ+1)+1:(2(ℓᵧ+1)+5)] .= Φ, δ, σ, δ_b, v_b  # write u with our variables
    u[2(ℓᵧ+1)+1:(2(ℓᵧ+1)+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end

function initial_conditions_spl(xᵢ, hierarchy::Hierarchy_spl{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    u = zeros(T, 2(ℓᵧ+1)+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    # Θ, Θᵖ, Φ, δ, σ, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ, Θᵖ, Φ, δ, v, δ_b, v_b = unpack_spl(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    #These get a 3/3 since massive neutrinos behave as massless at time of ICs
    f_ν = 1/(1 + 1/(7*(3/3)*3.046/8 *(4/11)^(4/3))) # we need to actually keep this
    α_c = par.α_c

    # metric and matter perturbations
    ℛ = 1.0  # set curvature perturbation to 1
    Φ = (4f_ν + 10) / (4f_ν + 15) * ℛ  # for a mode outside the horizon in radiation era
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) ) * Φ

    #trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓᵧ
        Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end

    δ = -α_c*Θ[0] # this is general enough to allow this to be any species
    δ_b = δ
    # v = -3k*Θ[1]
    v = α_c*k*Θ[1]
    v_b = -3k*Θ[1]

    u[2(ℓᵧ+1)+1:(2(ℓᵧ+1)+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end

#FIXME this is pretty old code that hasn't been tested in a while!
# TODO: this could be extended to any Newtonian gauge integrator if we specify the
# Bardeen potential Ψ and its derivative ψ′ for an integrator, or we saved them
function source_function(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute some quantities
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ′′(x)
    τₓ, τₓ′, τₓ′′ = ih.τ(x), ih.τ′(x), ih.τ′′(x)
    g̃ₓ, g̃ₓ′, g̃ₓ′′ = ih.g̃(x), ih.g̃′(x), ih.g̃′′(x)
    a = x2a(x)
    Θ, Θᵖ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ′, Θᵖ′, Φ′, δ′, v′, δ_b′, v_b′ = unpack(du, hierarchy)

    Ψ = -Φ - 12H₀² / k^2 / a^2 * (par.Ω_r * Θ[2]) #why am I doing this? - because H0 pulls out a factor of rho crit - just unit conversion
                                                                   #this introduces a factor of bg density I cancel using the integrated bg mnu density now

   Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * (par.Ω_r * (Θ′[2] - 2 * Θ[2]))

    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]

    term1 =  g̃ₓ * (Θ[0] + Ψ + Π/4) + exp(-τₓ) * (Ψ′ - Φ′)
    term2 = (-1/k) * (ℋₓ′ * g̃ₓ * v_b + ℋₓ * g̃ₓ′ * v_b + ℋₓ * g̃ₓ * v_b′)
    Π′′ = 2k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * Θ[1] + Θ′[1]) + (3/10) * (τₓ′′ * Π + τₓ′ * Π′) -
        3k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * (Θ[3] + Θᵖ[1] + Θᵖ[3]) + (Θ′[3] + Θᵖ′[1] + Θᵖ′[3]))
    term3 = (3/(4k^2)) * (
        (ℋₓ′^2 + ℋₓ * ℋₓ′′) * g̃ₓ * Π + 3 * ℋₓ * ℋₓ′ * (g̃ₓ′ * Π + g̃ₓ * Π′) +
        ℋₓ^2 * (g̃ₓ′′ * Π + 2g̃ₓ′ * Π′ + g̃ₓ * Π′′))
    return term1 + term2 + term3
end
