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

# wrapper for hierarchy that includes ctime conversion (and enforces type stability on η2x)
struct ConformalHierarchy{T<:Real,  H <: Hierarchy{T}, IT <: AbstractInterpolation{T}}
    hierarchy::H
    η2x::IT
end

function boltsolve(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6,abstol=1e-6) where T
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = initial_conditions(xᵢ, hierarchy)
    prob = ODEProblem{true}(hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,abstol=abstol,
                # saveat=hierarchy.bg.x_grid, 
                dense=true,
                )
    return sol
end

function boltsolve_conformal(confhierarchy::ConformalHierarchy{T},#FIXME we do't need this? {Hierarchy{T},AbstractInterpolation{T}},
                         ode_alg=KenCarp4(); reltol=1e-6,abstol=1e-6) where T
    hierarchy = confhierarchy.hierarchy
    xᵢ = hierarchy.bg.x_grid[1]#confhierarchy.η2x( hierarchy.bg.η(hierarchy.bg.x_grid[1]) )#η[1] ) #to be consistent
    u₀ = initial_conditions(xᵢ, hierarchy)
    # Mpcfac = hierarchy.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(hierarchy_conformal!, u₀, 
                            # (max(hierarchy.bg.η[1]*Mpcfac,hierarchy.bg.η(hierarchy.bg.x_grid[1])*Mpcfac), 
                            # min(hierarchy.bg.η[end]*Mpcfac,hierarchy.bg.η(hierarchy.bg.x_grid[end])*Mpcfac)),
                            (max(hierarchy.bg.η[1],hierarchy.bg.η(hierarchy.bg.x_grid[1])), 
                            min(hierarchy.bg.η[end],hierarchy.bg.η(hierarchy.bg.x_grid[end]))),
                            confhierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,abstol=abstol,
                # saveat=hierarchy.bg.η, 
                dense=true
                )
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

function hierarchy_conformal!(du, u, confhierarchy::ConformalHierarchy{T}, η) where T
    hierarchy = confhierarchy.hierarchy
    Mpcfac = hierarchy.bg.H₀*299792.458/100.
    x = confhierarchy.η2x(η)#  / Mpcfac )
    ℋ = hierarchy.bg.ℋ(x)
    hierarchy!(du, u, hierarchy, x)
    du .*= ℋ / Mpcfac  # account for dx/dη
    return nothing
end

# BasicNewtonian comes from Callin+06 and the Dodelson textbook (dispatches on hierarchy.integrator)
function hierarchy!(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    # q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    Mpcfac = hierarchy.bg.H₀*299792.458/100.
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x)/Mpcfac, ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)


    ℓ_ν = hierarchy.ℓ_ν
    ℓ_mν =  hierarchy.ℓ_mν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6

    # If using RSA, need to update the mono/quadrupole since they feed into metric perts
    τc = 1/(-τₓ′*ℋₓ)
    rsa_on = false#(k*ηₓ > 240) &&  (τc/ηₓ>100)
    if rsa_on
        Θ[0] = Φ - ℋₓ/k *τₓ′ * v_b
        Θ[2] = 0
        𝒩[0] = Φ
        𝒩[2] = 0
    end
    
    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩[2]#add rel quadrupole
                                  + σℳ / bg.ρ_crit /4
                                  ) 

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0] #add rel monopole on this line
        + a^(-2) * ρℳ / bg.ρ_crit
        ) 


    # if ((x<=-19.99 || x>=-0.01) &&  ~(typeof(Φ′) <: ForwardDiff.Dual))
    #     println("x = ", x)
    #     println("Φ′ = ", Φ′)
    #     println("𝒩[0] = ", 𝒩[0])
    #     println("Θ[0] = ", Θ[0])
    #     println("𝒩[2] = ", 𝒩[2])
    #     println("Ψ = ", Ψ)
    #     println("Ψ components: Θ₂ = $(Θ[2]), 𝒩₂ = $(𝒩[2]), σℳ = $(σℳ)")
    # end

    # RSA needs to come on first for baryons
    if rsa_on
        Θ[1] = ℋₓ/k * (  -2Φ′ + τₓ′*( Φ - csb²*δ_b  )
                         + ℋₓ/k*( τₓ′′ - τₓ′ )*v_b  )
        𝒩[1] = -2ℋₓ/k *Φ′
    end

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)

    # neutrinos (massive, MB 57)
    # for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
    for i_q in 0:nq-1
        q = xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        ϵ = √(q^2 + (a*m_ν)^2)
        df0 = dlnf0dlnq(q,par)
        #need these factors of 4 on Φ, Ψ terms due to MB pert defn
        ℳ′[0* nq+i_q] = - k / ℋₓ *  q/ϵ * ℳ[1* nq+i_q]  + Φ′ * df0
        ℳ′[1* nq+i_q] = k / (3ℋₓ) * ( q/ϵ * (ℳ[0* nq+i_q] - 2ℳ[2* nq+i_q])  - ϵ/q * Ψ  * df0)
        for ℓ in 2:(ℓ_mν-1)
            ℳ′[ℓ* nq+i_q] =  k / ℋₓ * q / ((2ℓ+1)*ϵ) * ( ℓ*ℳ[(ℓ-1)* nq+i_q] - (ℓ+1)*ℳ[(ℓ+1)* nq+i_q] )
        end
        ℳ′[ℓ_mν* nq+i_q] =  q / ϵ * k / ℋₓ * ℳ[(ℓ_mν-1)* nq+i_q] - (ℓ_mν+1)/(ℋₓ *ηₓ) *ℳ[(ℓ_mν)* nq+i_q] #MB (58) similar to rel case but w/ q/ϵ
    end

    # RSA equations (implementation of CLASS default switches)
    # rsa_on = (k*ηₓ > 240) &&  (-τₓ′*ηₓ*ℋₓ<100)
    #*sqrt(H₀²)< 1) #is this ℋ or H0?
    if rsa_on
        #photons
        # Θ[0] = Φ - ℋₓ/k *τₓ′ * v_b
        # Θ[1] = ℋₓ/k * (  -2Φ′ + τₓ′*( Φ - csb²*δ_b  )
        #                  + ℋₓ/k*( τₓ′′ - τₓ′ )*v_b  )
        # Θ[2] = 0

        #massless neutrinos
        # 𝒩[0] = Φ
        # 𝒩[1] = -2ℋₓ/k *Φ′
        # 𝒩[2] = 0

        #set polarization to zero
        Θᵖ[0] = 0
        Θᵖ[1] = 0
        Θᵖ[2] = 0

        # manual zeroing to avoid saving garbage
        𝒩′[:] = zeros(ℓ_ν+1)
        Θ′[:] = zeros(ℓᵧ+1)
        Θᵖ′[:] = zeros(ℓᵧ+1)

    else
        #do usual hierarchy
        # relativistic neutrinos (massless)
        𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
        𝒩′[1] = k/(3ℋₓ) * 𝒩[0] - 2*k/(3ℋₓ) *𝒩[2] + k/(3ℋₓ) *Ψ
        for ℓ in 2:(ℓ_ν-1)
            𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) * ( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1] )
        end
        #truncation (MB/Dodelson)
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

    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end

# BasicNewtonian Integrator (dispatches on hierarchy.integrator)
function initial_conditions(xᵢ, hierarchy::Hierarchy{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih, hierarchy.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    # q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    ℓ_ν = hierarchy.ℓ_ν
    ℓ_mν =  hierarchy.ℓ_mν
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)
    Mpcfac = hierarchy.bg.H₀*299792.458/100.
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ)/Mpcfac, ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    #These get a 3/3 since massive neutrinos behave as massless at time of ICs
    f_ν = 1/(1 + 1/(7*(3/3)*par.N_ν/8 *(4/11)^(4/3)))

    # metric and matter perturbations
    Φ = 1.0 #-0.0008000688458547067*par.h #(1 + 2/5 * f_ν) / (3/2 + 2/5 * f_ν) / par.h #1.0
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) )

    #trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k * ηₓ) / 3 # this was for clarity but wastes  (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1] #This is not in MB (and is not strictly consistent)...idea is should be effectively zero since tau' huge early on...
    
    #^Numerically this is irrelevant at superhorizon (it is 10⁻⁸𝒩[2]), but allegedly helps with numerics...
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

    
    # neutrino hierarchy
    # we need xᵢ to be before neutrinos decouple, as always
    𝒩[0] = Θ[0]
    𝒩[1] = Θ[1]
    𝒩[2] = - (k^2 *ηₓ^2)/15 * 1 / (1 + 2/5 *f_ν) * Φ  / 2 #MB
    #FIXME^put the C here for consistency
    # 𝒩[2] = - (k^2 *aᵢ²*Φ) / (12H₀² * Ω_ν) * 1 / (1 + 5/(2*f_ν)) #Callin06
    #These are the same to 3 decimal places ...about the expected error on conformal time spline
    for ℓ in 3:ℓ_ν
        𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #standard truncation
    end

    #massive neutrino hierarchy
    #It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    # for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
    for i_q in 0:nq-1
        q = xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        ϵ = √(q^2 + (aᵢ*par.Σm_ν)^2)
        df0 = dlnf0dlnq(q,par)
        ℳ[0* nq+i_q] = -𝒩[0]  *df0
        ℳ[1* nq+i_q] = -ϵ/q * 𝒩[1] *df0
        ℳ[2* nq+i_q] = -𝒩[2]  *df0  #drop quadratic+ terms in (ma/q) as in MB
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
    ρ0ℳ = bg.ρ₀ℳ(x) #get current value of massive neutrino backround density from spline
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    Ω_ν =  7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    # q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, Φ′, δ′, v′, δ_b′, v_b′ = unpack(du, hierarchy)

    # recalulate these since we didn't save them (Callin eqns 39-42)
    #^Also have just copied from before, but should save these maybe?
    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    _, σℳ′ = ρ_σ(ℳ′[0:nq-1], ℳ′[2*nq:3*nq-1], bg, a, par)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩[2]#add rel quadrupole
                                  + σℳ / bg.ρ_crit /4
                                  )

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

# The polarization source function (jms 6/6/22 UNTESTED!) SZ eqn 12d (in our units)
function source_function_P(du, u, hierarchy::Hierarchy{T, BasicNewtonian}, x) where T
    k, ℓᵧ, par, bg, ih,nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ′′(x)
    τₓ, τₓ′, τₓ′′ = ih.τ(x), ih.τ′(x), ih.τ′′(x)
    g̃ₓ, g̃ₓ′, g̃ₓ′′ = ih.g̃(x), ih.g̃′(x), ih.g̃′′(x)
    a = x2a(x)
    ρ0ℳ = bg.ρ₀ℳ(x) #get current value of massive neutrino backround density from spline
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    Ω_ν =  7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    # q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ are mutable views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, Φ′, δ′, v′, δ_b′, v_b′ = unpack(du, hierarchy)


    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    return (3/(4k^2)) * g̃ₓ * Π 
end

function ρ_σ(ℳ0,ℳ2,bg,a,par::AbstractCosmoParams) #a mess
    #Do q integrals to get the massive neutrino metric perturbations
    #MB eqn (55)
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    #^Replace this with bg.ρ_crit? I think it is using an imported function ρ_crit
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)

    #FIXME: avoid repeating code? and maybe put general integrals in utils?
    m = par.Σm_ν
    nq = length(ℳ0) #assume we got this right
    ϵx(x, am) = √(xq2q(x,logqmin,logqmax)^2 + (am)^2)
    Iρ(x) = xq2q(x,logqmin,logqmax)^2  * ϵx(x, a*m) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    Iσ(x) = xq2q(x,logqmin,logqmax)^2  * (xq2q(x,logqmin,logqmax)^2 /ϵx(x, a*m)) * f0(xq2q(x,logqmin,logqmax),par) / dxdq(xq2q(x,logqmin,logqmax),logqmin,logqmax)
    xq,wq = bg.quad_pts,bg.quad_wts
    # ρ = 4π*sum(Iρ.(xq).*ℳ0.*wq)
    # σ = 4π*sum(Iσ.(xq).*ℳ2.*wq)
    ρ,σ=0.0,0.0
    for i in 1:length(xq)
        ρ+=Iρ(xq[i])*ℳ0[i]*wq[i] #confusingly, since ℳ0 is a view it starts at 1, not actually offset array...
        σ+=Iσ(xq[i])*ℳ2[i]*wq[i]
    end

    # #a-dependence has been moved into Einstein eqns, as have consts in σ
    return 4π*ρ,4π*σ
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

function rsa_perts!(u, hierarchy::Hierarchy{T},x) where T
    #redundant code for what we need to compute RSA perts in place in u
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)
    ℓ_ν = hierarchy.ℓ_ν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    # Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6

    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩[2]
                                  + σℳ / bg.ρ_crit /4
                                  )
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
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

    #set polarization to zero
    Θᵖ[0] = 0
    Θᵖ[1] = 0
    Θᵖ[2] = 0

    u[1] = Θ[0]
    u[2] = Θ[1]
    u[3] = Θ[2]

    u[(ℓᵧ+1)+1] = Θᵖ[0]
    u[(ℓᵧ+1)+2] = Θᵖ[1]
    u[(ℓᵧ+1)+3] = Θᵖ[2]

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

function boltsolve_rsa(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6) where T
    #call solve as usual first
    perturb = boltsolve(hierarchy, reltol=reltol)
    x_grid = hierarchy.bg.x_grid
    pertlen = 2(hierarchy.ℓᵧ+1)+(hierarchy.ℓ_ν+1)+(hierarchy.ℓ_mν+1)*hierarchy.nq+5
    results=zeros(pertlen,length(x_grid))
    for i in 1:length(x_grid) results[:,i] = perturb(x_grid[i]) end
    #replace the late-time perts with RSA approx (assuming we don't change rsa switch)
    # this_rsa_switch = x_grid[argmin(abs.(hierarchy.k .* hierarchy.bg.η.(x_grid) .- 45))]
    #⬇check if we are always outside horizon, if so, then never turn on rsa
    xrsa_hor =  sum((@. hierarchy.k*hierarchy.bg.η .> 240)) > 0  ? minimum(x_grid[(@. hierarchy.k*hierarchy.bg.η .> 240)]) : x_grid[end]
    xrsa_od = minimum(x_grid[(@. -hierarchy.ih.τ′*hierarchy.bg.η*hierarchy.bg.ℋ .<100)])
    this_rsa_switch = max(xrsa_hor,xrsa_od)
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