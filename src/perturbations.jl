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



function boltsolve(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6, abstol=1e-6) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, abstol=abstol,
                saveat=hierarchy.bg.x_grid, dense=false,
                )
    return sol
end

function rsa_perts!(u, hierarchy::Hierarchy{T},x) where T
    #redundant code for what we need to compute RSA perts in place in u
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    Ω_r, Ω_b, Ω_c, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)
    ℓ_ν = hierarchy.ℓ_ν
    (Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b) = unpack(u, hierarchy)

    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩[2]
                                  + σℳ / bg.ρ_crit /4
                                  )
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0]
        + a^(-2) * ρℳ / bg.ρ_crit
        )

    #fixed RSA
    Θ[0] = Φ - ℋₓ/k *τₓ′ * v_b
    Θ[1] = ℋₓ/k * (  -2Φ′ + τₓ′*( Φ - csb²*δ_b  )
                     + ℋₓ/k*( τₓ′′ - τₓ′ )*v_b  )
    Θ[2] = 0
    #massless neutrinos
    𝒩[0] = Φ
    𝒩[1] = -2ℋₓ/k *Φ′
    𝒩[2] = 0

    u[1] = Θ[0]
    u[2] = Θ[1]
    u[3] = Θ[2]

    u[2(ℓᵧ+1)+1] = 𝒩[0]
    u[2(ℓᵧ+1)+2] = 𝒩[1]
    u[2(ℓᵧ+1)+3] = 𝒩[2]

    #zero the rest to avoid future confusion
    for ℓ in 3:(ℓᵧ)
        u[ℓ] = 0
        u[(ℓᵧ+1)+ℓ] = 0
    end
    for ℓ in 3:(ℓ_ν) u[2(ℓᵧ+1)+ℓ] = 0 end
    return nothing
end

function boltsolve_rsa(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6, abstol=1e-6) where T
    #call solve as usual first
    perturb = boltsolve(hierarchy, reltol=reltol, abstol=abstol)
    x_grid = hierarchy.bg.x_grid
    pertlen = 2(hierarchy.ℓᵧ+1)+(hierarchy.ℓ_ν+1)+(hierarchy.ℓ_mν+1)*hierarchy.nq+5
    results=zeros(pertlen,length(x_grid))
    for i in 1:length(x_grid) results[:,i] = perturb(x_grid[i]) end
    #replace the late-time perts with RSA approx (assuming we don't change rsa switch)
    #this_rsa_switch = x_grid[argmin(abs.(hierarchy.k .* hierarchy.bg.η.(x_grid) .- 45))]

    xrsa_hor = findfirst(>(240), @. hierarchy.k * hierarchy.bg.η)
    xrsa_od = findfirst(>(100), @. -hierarchy.ih.τ′*hierarchy.bg.ℋ/hierarchy.bg.η)
    xrsa_hor = isnothing(xrsa_hor) ? length(x_grid) : xrsa_hor
    xrsa_od = isnothing(xrsa_hor) ? length(x_grid) : xrsa_od

    this_rsa_switch = x_grid[max(xrsa_hor,xrsa_od)]
    x_grid_rsa = x_grid[x_grid.>this_rsa_switch]
    results_rsa = results[:,x_grid.>this_rsa_switch]
    #(re)-compute the RSA perts so we can write them to the output vector
    for i in 1:length(x_grid_rsa)
        rsa_perts!(view(results_rsa,:,i),hierarchy,x_grid_rsa[i]) #to mutate need to use view...
    end
    results[:,x_grid.>this_rsa_switch] = results_rsa
    sol = results
    return sol
end

# basic Newtonian gauge: establish the order of perturbative variables in the ODE solve
function unpack(u, ::Hierarchy{<:Any,BasicNewtonian})
    Θ = Origin(0)(u.Θ)
    Θᵖ = Origin(0)(u.Θᵖ)
    𝒩 = Origin(0)(u.𝒩)
    ℳ = Origin(0,1)(u.ℳ)
    @unpack (Φ, δ, v, δ_b, v_b) = u
    return (Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b)
end

function ρ_σ(ℳ0, ℳ2, bg, a, par::AbstractCosmoParams)
    # Do q integrals to get the massive neutrino metric perturbations
    # MB eqn (55)
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 * bg.ρ_crit * par.Ω_r)^(1/4)
    logqmin, logqmax = log10(Tν/30), log10(Tν*30)
    # FIXME: avoid repeating code? and maybe put general integrals in utils?
    m = par.Σm_ν
    ϵx(x, am) = √(xq2q(x,logqmin,logqmax)^2 + (am)^2)
    Iρ(x) = xq2q(x,logqmin,logqmax)^2 * ϵx(x, a*m) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    Iσ(x) = xq2q(x,logqmin,logqmax)^2 * (xq2q(x,logqmin,logqmax)^2 /ϵx(x, a*m)) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    ρ = σ = zero(Tν)
    for qᵢ in 1:length(bg.quad_pts)
        ρ += Iρ(bg.quad_pts[qᵢ]) * ℳ0[qᵢ] * bg.quad_wts[qᵢ]
        σ += Iσ(bg.quad_pts[qᵢ]) * ℳ2[qᵢ] * bg.quad_wts[qᵢ]
    end
    # #a-dependence has been moved into Einstein eqns, as have consts in σ
    return 4π*ρ, 4π*σ
end

#need a separate function for θ (really(ρ̄+P̄)θ) for plin gauge change
function θ(ℳ1,bg,a,par::AbstractCosmoParams) #a mess
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *bg.ρ_crit *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    m = par.Σm_ν
    nq = length(ℳ1) #assume we got this right
    Iθ(x) = xq2q(x,logqmin,logqmax)^3  * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    xq,wq = bg.quad_pts,bg.quad_wts
    θ = 4π*sum(Iθ.(xq).*ℳ1.*wq)
    #Note that this still needs to be multiplied with ka^-4 prefactor
    return θ
end

# BasicNewtonian comes from Callin+06 and the Dodelson textbook (dispatches on hierarchy.integrator)
function hierarchy!(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    @unpack (k, ℓᵧ, ℓ_ν, ℓ_mν, par, bg, ih, nq) = hierarchy
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 * bg.ρ_crit * par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    Ω_r, Ω_b, Ω_c, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)

    (Θ,  Θᵖ,  𝒩,  ℳ, Φ, δ, v, δ_b, v_b) = unpack(u,  hierarchy)
    (Θ′, Θᵖ′, 𝒩′, ℳ′)                   = unpack(du, hierarchy)


    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ = @views ρ_σ(ℳ[0,:], ℳ[2,:], bg, a, par) # monopole (energy density, 00 part), quadrupole (shear stress, ij part)
    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩[2]#add rel quadrupole
                                  + σℳ / bg.ρ_crit /4
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0] #add rel monopole on this line
        + a^(-2) * ρℳ / bg.ρ_crit
        )

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)

    # neutrinos (massive, MB 57)
    # TODO: it might be possible to transpose this loop and get it hardware vectorized
    for qᵢ in 1:length(bg.quad_pts)
        q = xq2q(bg.quad_pts[qᵢ], logqmin, logqmax)
        ϵ = √(q^2 + (a*m_ν)^2)
        df0 = dlnf0dlnq(q, par)
        # need these factors of 4 on Φ, Ψ terms due to MB pert defn
        ℳ′[0,qᵢ] = - k / ℋₓ *  q/ϵ * ℳ[1,qᵢ]  + Φ′ * df0
        ℳ′[1,qᵢ] = k / (3ℋₓ) * ( q/ϵ * (ℳ[0,qᵢ] - 2ℳ[2,qᵢ])  - ϵ/q * Ψ  * df0)
        for ℓ in 2:(ℓ_mν-1)
            ℳ′[ℓ,qᵢ] =  k / ℋₓ * q / ((2ℓ+1)*ϵ) * ( ℓ*ℳ[ℓ-1,qᵢ] - (ℓ+1)*ℳ[ℓ+1,qᵢ] )
        end
        ℳ′[ℓ_mν,qᵢ] =  q / ϵ * k / ℋₓ * ℳ[ℓ_mν-1,qᵢ] - (ℓ_mν+1)/(ℋₓ *ηₓ) *ℳ[ℓ_mν,qᵢ] #MB (58) similar to rel case but w/ q/ϵ
    end

    # RSA equations (implementation of CLASS default switches)
    rsa_on = (k*ηₓ > 240) &&  (-τₓ′*ℋₓ / ηₓ > 100)
    if rsa_on
        # photons
        Θ[0] = Φ - ℋₓ/k *τₓ′ * v_b
        Θ[1] = -2Φ′/k + (k^-2)*( τₓ′′ * v_b + τₓ′ * (ℋₓ*v_b - csb² *δ_b/k + k*Φ) )
        Θ[1] = ℋₓ/k * (  -2Φ′ + τₓ′*( Φ - csb²*δ_b  )
                         + ℋₓ/k*( τₓ′′ - τₓ′ )*v_b  )
        Θ[2] = 0
        # massless neutrinos
        𝒩[0] = Φ
        𝒩[1] = -2ℋₓ/k *Φ′
        𝒩[2] = 0

        # manual zeroing to avoid saving garbage
        𝒩′ .= 0
        Θ′ .= 0
        Θᵖ′ .= 0

    else
        #do usual hierarchy
        # relativistic neutrinos (massless)
        𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
        𝒩′[1] = k/(3ℋₓ) * 𝒩[0] - 2*k/(3ℋₓ) *𝒩[2] + k/(3ℋₓ) *Ψ
        for ℓ in 2:(ℓ_ν-1)
            𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) * ( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1] )
        end
        #truncation (same between MB and Callin06/Dodelson)
        𝒩′[ℓ_ν] =  k / ℋₓ  * 𝒩[ℓ_ν-1] - (ℓ_ν+1)/(ℋₓ *ηₓ) *𝒩[ℓ_ν]


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

    end
    #END RSA

    (Φ=Φ′; δ=δ′; v=v′; δ_b=δ_b′; v_b=v_b′)
    @pack! du = (Φ, δ, v, δ_b, v_b)
    
    return nothing
end

# BasicNewtonian Integrator (dispatches on hierarchy.integrator)
function initial_conditions(xᵢ, hierarchy::Hierarchy{T,BasicNewtonian}) where {T}

    @unpack (k, ℓᵧ, ℓ_ν, ℓ_mν, par, bg, ih, nq) = hierarchy

    Tν = (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/π^2 * bg.ρ_crit * par.Ω_r)^(1/4)
    (logqmin, logqmax) = log10(Tν/30), log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts, logqmin, logqmax)

    u = ComponentVector{T}(Φ=0, δ=0, v=0, δ_b=0, v_b=0, Θ=zeros(ℓᵧ+1), Θᵖ=zeros(ℓᵧ+1), 𝒩=zeros(ℓ_ν+1), ℳ=zeros(ℓ_mν+1,nq))

    ℋₓ, ηₓ, τₓ′ = bg.ℋ(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ)
    (Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b) = unpack(u, hierarchy)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    # These get a 3/3 since massive neutrinos behave as massless at time of ICs
    Ω_ν = (7/8) * (3/3) * par.N_ν * (4/11)^(4/3) * par.Ω_r
    f_ν = 1 / (1 + 1 / ((7/8) * (3/3) * par.N_ν * (4/11)^(4/3)))
    # ρ0ℳ = bg.ρ₀ℳ(xᵢ)

    # metric and matter perturbations
    ℛ = 1.0  # set curvature perturbation to 1
    Φ = (4f_ν + 10) / (4f_ν + 15) * ℛ  # for a mode outside the horizon in radiation era
    C = -((15 + 4f_ν) / (20 + 8f_ν)) * Φ

    # trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓᵧ
        Θ[ℓ]  = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end

    δ = 3/4 * 4Θ[0] # the 4 converts δγ_MB -> Dodelson convention
    δ_b = δ
    # we have that Θc = Θb = Θγ = Θν, but need to convert Θ = - k v (i absorbed in v)
    v = -3k * Θ[1]
    v_b = v

    # neutrino hierarchy
    # we need xᵢ to be before neutrinos decouple, as always
    𝒩[0] = Θ[0]
    𝒩[1] = Θ[1]
    𝒩[2] = -(k^2 * ηₓ^2)/15 * 1 / (1 + 2/5 * f_ν) * Φ / 2 #MB
    #FIXME^put the C here for consistency
    # println("MB nu quad: ", - (k^2 *ηₓ^2)/30 * 1 / (1 + 2/(5) *f_ν) * Φ)
    # println("Callin nu quad ", - (k^2 *aᵢ²*Φ) / (12H₀² * Ω_ν) * 1 / (1 + 5/(2*f_ν)))
    # 𝒩[2] = - (k^2 *aᵢ²*Φ) / (12H₀² * Ω_ν) * 1 / (1 + 5/(2*f_ν)) #Callin06
    #These are the same to 3 decimal places ...about the expected error on conformal time spline
    for ℓ in 3:ℓ_ν
        𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #standard truncation
    end

    # massive neutrino hierarchy
    # It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    for (qᵢ, q) in enumerate(q_pts)
        ϵ = √(q^2 + (aᵢ*par.Σm_ν)^2)
        df0 = dlnf0dlnq(q, par)
        ℳ[0,qᵢ] = -𝒩[0] * df0
        ℳ[1,qᵢ] = -ϵ/q * 𝒩[1] *df0
        ℳ[2,qᵢ] = -𝒩[2] * df0  #drop quadratic+ terms in (ma/q) as in MB
        for ℓ in 3:ℓ_mν #same scheme for higher-ell as for relativistic
            ℳ[ℓ,qᵢ] = q / ϵ * k/((2ℓ+1)ℋₓ) * ℳ[ℓ-1,qᵢ] #approximation equivalent to MB, but add q/ϵ - leaving as 0 makes no big difference
        end
    end

    @pack! u = (Φ, δ, v, δ_b, v_b)

    return u

end

#FIXME this is pretty old code that hasn't been tested in a while!
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
    ρ0ℳ = bg.ρ₀ℳ(x) #get current value of massive neutrino backround density from spline
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 * bg.ρ_crit * par.Ω_r)^(1/4)
    Ω_ν =  7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    (Θ,  Θᵖ,  𝒩,  ℳ,  Φ,  δ,  v,  δ_b,  v_b ) = unpack(u,  hierarchy)
    (Θ′, Θᵖ′, 𝒩′, ℳ′, Φ′, δ′, v′, δ_b′, v_b′) = unpack(du, hierarchy)

    # recalulate these since we didn't save them (Callin eqns 39-42)
    #^Also have just copied from before, but should save these maybe?
    _, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    _, σℳ′ = ρ_σ(ℳ′[0:nq-1], ℳ′[2*nq:3*nq-1], bg, a, par)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (par.Ω_r * Θ[2]
                                  + Ω_ν * 𝒩[2] #add rel quadrupole
                                  + σℳ / bg.ρ_crit) #why am I doing this? - because H0 pulls out a factor of rho crit - just unit conversion
                                                                   #this introduces a factor of bg density I cancel using the integrated bg mnu density now

   Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * (par.Ω_r * (Θ′[2] - 2 * Θ[2])
                                   + Ω_ν * (𝒩′[2] - 2 * 𝒩[2])
                                   + (σℳ′ - 2 * σℳ) / bg.ρ_crit /4 )

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
