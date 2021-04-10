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
    ℓ_ν::Int
    ℓ_mν::Int
    nq::Int
end

Hierarchy(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real, ℓᵧ=8, ℓ_ν=8, ℓ_mν=10, nq=15) = Hierarchy(integrator, par, bg, ih, k, ℓᵧ, ℓ_ν,ℓ_mν, nq)

function boltsolve(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=hierarchy.bg.x_grid, dense=false)
    return sol
end

# basic Newtonian gauge: establish the order of perturbative variables in the ODE solve
function unpack(u, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    ℓᵧ = hierarchy.ℓᵧ
    ℓ_ν =  hierarchy.ℓ_ν
    ℓ_mν = hierarchy.ℓ_mν #should be smaller than others
    nq = hierarchy.nq
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    𝒩 = OffsetVector(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1)) , 0:ℓ_ν)  # indexed 0 through ℓ_ν
    ℳ = OffsetVector(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq )) , 0:(ℓ_mν+1)*nq -1)  # indexed 0 through ℓ_mν
    Φ, δ, v, δ_b, v_b = view(u, ((2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5)) #getting a little messy...
    return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
end

function ρ_σ(ℳ0,ℳ2,bg,a,par::AbstractCosmoParams) #a mess
    #Do q integrals to get the massive neutrino metric perturbations
    #MB eqn (55)
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)#1e-6,1e-1

    #FIXME: avoid repeating code? and maybe put general integrals in utils?
    m = par.Σm_ν
    nq = length(ℳ0) #assume we got this right
    ϵx(x, am) = √(xq2q(x,logqmin,logqmax)^2 + (am)^2)
    Iρ(x) = xq2q(x,logqmin,logqmax)^2  * ϵx(x, a*m) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    Iσ(x) = xq2q(x,logqmin,logqmax)^2  * (xq2q(x,logqmin,logqmax)^2 /ϵx(x, a*m)) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)

    xq,wq = bg.quad_pts,bg.quad_wts
    ρ = 4π*sum(Iρ.(xq).*ℳ0.*wq)
    σ = 4π*sum(Iσ.(xq).*ℳ2.*wq)
    # #a-dependence has been moved into Einstein eqns, as have consts in σ
    return ρ,σ
end

# BasicNewtonian comes from Callin+06 and the Dodelson textbook (dispatches on hierarchy.integrator)
function hierarchy!(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    ℓ_ν = hierarchy.ℓ_ν
    ℓ_mν =  hierarchy.ℓ_mν
    norm𝒩′ = 1.0 /(Ω_ν * bg.ρ_crit / 2)# par.N_ν) #Normalization to match 𝒩 after integrating, par.N_ν->2
    norm𝒩 = norm𝒩′/ 4.0
    #^Here we remove the 4 in denom b/c it has moved to the Einstein eqns.

    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6

    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                  + Ω_ν * 𝒩[2] #add rel quadrupole
                                  + σℳ / bg.ρ_crit/ norm𝒩′) #add mnu integrated quadrupole

    # println("x= ",x, " so a = ", exp(x))
    # println("Size of terms in i neq j eqn. Ω_ν: ", Ω_ν * 𝒩[2], " and σℳ ", σℳ / bg.ρ_crit / norm𝒩′)

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0] #add rel monopole on this line
        + a^(-2) * ρℳ  / bg.ρ_crit / norm𝒩′) #add mnu integrated monopole

    # println("Size of terms in 00 eqn. Ω_ν: ", 4Ω_ν * a^(-2) * 𝒩[0], " and ρℳ ", 4 * a^(-2) * ρℳ  / bg.ρ_crit / norm𝒩′)

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)

    # relativistic neutrinos (massless)
    𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
    𝒩′[1] = k/(3ℋₓ) * 𝒩[0] - 2*k/(3ℋₓ) *𝒩[2] + k/(3ℋₓ) *Ψ
    for ℓ in 2:(ℓ_ν-1) #ℓ_ν same as ℓᵧ for massless nu for now
        𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) * ( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1] )
    end
    #truncation (same between MB and Callin06)
    𝒩′[ℓ_ν] =  k / ℋₓ  * 𝒩[ℓ_ν-1] - (ℓ_ν+1)/(ℋₓ *ηₓ) *𝒩[ℓ_ν]

    #WIP: nonrelativistic nu
    # neutrinos (massive, MB 57)
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
        ϵ = √(q^2 + (a*m_ν)^2)
        #dlnf0dlnq = bg.df0(log10(q)) * norm𝒩
        df0 = dlnf0dlnq(q,par) * norm𝒩
        #need these factors of 4 on Φ, Ψ terms due to MB pert defn
        ℳ′[0* nq+i_q] = - k / ℋₓ *  q/ϵ * ℳ[1* nq+i_q] + Φ′ * df0
        ℳ′[1* nq+i_q] = k / (3ℋₓ) * (( q/ϵ * (ℳ[0* nq+i_q] - 2ℳ[2* nq+i_q])) - ϵ/q * Ψ  * df0)
        for ℓ in 2:(ℓ_mν-1)
            ℳ′[ℓ* nq+i_q] =  k / ℋₓ * q / ((2ℓ+1)*ϵ) * ( ℓ*ℳ[(ℓ-1)* nq+i_q] - (ℓ+1)*ℳ[(ℓ+1)* nq+i_q] )
        end
        ℳ′[ℓ_mν* nq+i_q] =  q / ϵ * k / ℋₓ * ℳ[(ℓ_mν-1)* nq+i_q] - (ℓ_mν+1)/(ℋₓ *ηₓ) *ℳ[(ℓ_mν)* nq+i_q] #MB (58) similar to rel case but w/ q/ϵ
    end

    #check monopole, dipole, quadrupole
    # ρℳ′, σℳ′  =  ρ_σ(ℳ′[0:nq-1], ℳ′[2*nq:3*nq-1], bg, a, par)
    # println("Size of 𝒩0` : ",𝒩′[0] , " and ρℳ` ",  ρℳ′)
    # println("Size of 𝒩2` : ",𝒩′[2] , " and σℳ` ",  σℳ′)
    # θℳ′, _  =  ρ_σ(ℳ′[nq:2*nq-1], zeros(nq), bg, a, par) #approximate ϵ=q
    # println("Size of 𝒩1` : ",𝒩′[1] , " and θℳ` ",  θℳ′)
    # maxℳ′, _  =  ρ_σ(ℳ′[(ℓ_mν-1)*nq:ℓ_mν*nq-1], zeros(nq), bg, a, par) #not sure if kosher
    # println("Size of max1` : ",𝒩′[ℓ_ν] , " and maxℳ` ",  maxℳ′)
    #
    # #check sizes of individual terms
    # println("Φ′ term - massless: ", -Φ′)
    # df0test = [dlnf0dlnq(q,par) for q in q_pts]
    # println("Φ′ term - massive: ", Φ′ * ρ_σ(df0test * norm𝒩, zeros(nq), bg, a, par)[1])
    # println("Ψ term - massless: ",k/(3ℋₓ) *Ψ)
    # println("Ψ term - massive: ",k/(3ℋₓ)* Ψ *ρ_σ(- sqrt.(ones(nq) .+ (a*m_ν ./ q_pts).^2)  .* df0test * norm𝒩, zeros(nq), bg, a, par)[1])


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

    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end

# BasicNewtonian Integrator (dispatches on hierarchy.integrator)
function initial_conditions(xᵢ, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih, hierarchy.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    ℓ_ν = hierarchy.ℓ_ν
    ℓ_mν =  hierarchy.ℓ_mν
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)

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
    # we need xᵢ to be before neutrinos decouple
    Ω_ν =  7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
    f_ν = 1/(1 + 1/(7*(2/3)*par.N_ν/8 *(4/11)^(4/3)))
    𝒩[0] = Θ[0]
    𝒩[1] = Θ[1]
    𝒩[2] = - (k^2 *aᵢ²*Φ) / (12H₀² * Ω_ν) * 1 / (1 + 5/(2*f_ν)) #Callin06 (71)
    for ℓ in 3:ℓ_ν
        𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #approximation of Callin06 (72)
    end

    #massive neutrino hierarchy
    #It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    norm𝒩 = 1/(4Ω_ν * bg.ρ_crit / 2)#par.N_ν) #Normalization to match 𝒩 after integrating, par.N_ν->2
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
        ϵ = √(q^2 + (aᵢ*par.Σm_ν)^2)
        df0 = dlnf0dlnq(q,par) * norm𝒩
        ℳ[0* nq+i_q] = -𝒩[0]  *df0
        ℳ[1* nq+i_q] = -ϵ/q * 𝒩[1] *df0
        ℳ[2* nq+i_q] = -𝒩[2]  *df0 #drop quadratic+ terms in (ma/q) as in MB
        for ℓ in 3:ℓ_mν #same scheme for higher-ell as for relativistic
            ℳ[ℓ* nq+i_q] = q / ϵ * k/((2ℓ+1)ℋₓ) * ℳ[(ℓ-1)*nq+i_q] #approximation of Callin06 (72), but add q/ϵ - leaving as 0 makes no big difference
        end
    end

    u[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end

# TODO: this could be extended to any Newtonian gauge integrator if we specify the
# Bardeen potential Ψ and its derivative ψ′ for an integrator, or we saved them
function source_function(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute some quantities
    k, ℓᵧ, par, bg, ih,nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ′′(x)
    τₓ, τₓ′, τₓ′′ = ih.τ(x), ih.τ′(x), ih.τ′′(x)
    g̃ₓ, g̃ₓ′, g̃ₓ′′ = ih.g̃(x), ih.g̃′(x), ih.g̃′′(x)
    a = x2a(x)

    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, Φ′, δ′, v′, δ_b′, v_b′ = unpack(du, hierarchy)

    # recalulate these since we didn't save them (Callin eqns 39-42)
    #FIXME check the neutrino contributions to Ψ and Ψ′!
    #^Also have just copied from before, but should save these maybe?
    Ω_ν =  7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *Ω_r
    norm𝒩 = 1/(4Ω_ν * bg.ρ_crit / par.N_ν)
    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    _, σℳ′ = ρ_σ(ℳ′[0:nq-1], ℳ′[2*nq:3*nq-1], bg, a, par)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                  + Ω_ν * 𝒩[2] #add rel quadrupole
                                  + σℳ / bg.ρ_crit / norm𝒩) #add mnu integrated quadrupole

    Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * (par.Ω_r * (Θ′[2] - 2 * Θ[2])
                                    + Ω_ν * (𝒩′[2] - 2 * 𝒩[2])
                                    + (σℳ′ - 2 * σℳ) / bg.ρ_crit/ norm𝒩)
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
