# these types and functions integrate the Boltzmann hierarchy through time

#for now just for photons, swap ℓᵧ for Nᵧ the number of IE conformal time integration points
#at the moment neutrinos are still hierarchy, will eventually replace neutrinos as well
struct IE{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
                 BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real,
				 IT<:AbstractInterpolation{T,1}}
    integrator::PI
    par::CP
    bg::BG
    ih::IH
    k::Tk
    sΘ2::IT
    sΠ::IT
    Nᵧ₁::Int #pre-entry
    Nᵧ₂::Int #recomb
    Nᵧ₃::Int #post-recomb
    ℓ_ν::Int
    ℓ_mν::Int
    nq::Int
end

# TODO for now I am copying for a separate neutrino integrator to simplify testing,
# but need to put these together (i.e. no more ℓ parameters, only splines).
# Massive neutrinos will have to be a collection of splines somehow...custom type?
# For now this only does massless neutrinos and evolves the rest via hierarchy
struct IEν{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
    BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real,
    IT<:AbstractInterpolation{T,1}}
    integrator::PI
    par::CP
    bg::BG
    ih::IH
    k::Tk
    s𝒩₀::IT
    s𝒩₂::IT
    N_ν::Int
    ℓ_γ::Int
    ℓ_mν::Int
    nq::Int
end

IE(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real,
    sΘ2::AbstractInterpolation,sΠ::AbstractInterpolation,
	#^FIXME: Are these right?? I dropped the {T,1} since T is not known here to get it to compile
    Nᵧ₁=10, Nᵧ₂=100, Nᵧ₃=50,
    ℓ_ν=8, ℓ_mν=10, nq=15
    ) = IE(integrator, par, bg, ih, k, sΘ2, sΠ, Nᵧ₁,Nᵧ₂,Nᵧ₃, ℓ_ν,ℓ_mν, nq)

IEν(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real,
    s𝒩₀::AbstractInterpolation,s𝒩₂::AbstractInterpolation,
    N_ν=400, ℓ_γ=8, ℓ_mν=10, nq=15
    ) = IEν(integrator, par, bg, ih, k, s𝒩₀, s𝒩₂, N_ν, ℓ_γ,ℓ_mν, nq)



function itersolve(Nₖ::Int,ie_0::IE{T};reltol=1e-6) where T
    x_grid = x_grid_ie(ie_0)
    Θ₂,Π =  zeros(length(x_grid)),zeros(length(x_grid)) #initialize to zero (for now)
    pertlen = 2(2+1)+(ie_0.ℓ_ν+1)+(ie_0.ℓ_mν+1)*ie_0.nq+5
    u_all = zeros(pertlen,length(x_grid))
    for k in 1:Nₖ
            Θ₂,Π,u_all = iterate(Θ₂,Π, ie_0.par, ie_0.bg, ie_0.ih, ie_0.k,  
                                    ie_0.Nᵧ₁,ie_0.Nᵧ₂,ie_0.Nᵧ₃, 
                                    x_grid, ie_0.ℓ_ν, ie_0.ℓ_mν, ie_0.nq, 
                                    reltol)
    end
    return u_all
end

function x_grid_ie(ie) 
    bg,ih,k = ie.bg,ie.ih,ie.k
    # Three phases: 
    # 1. Pre-horizon entry:
    xhor = bg.x_grid[argmin(abs.(k .* bg.η .- 1))] #horizon crossing ish
    x_ph_i, x_ph_f, n_ph = bg.x_grid[1], xhor, ie.Nᵧ₁ #10.
    dx_ph = (x_ph_f-x_ph_i)/(n_ph-1)
    x_ph = -20.:dx_ph:x_ph_f
    x_ph = -20. .+ dx_ph*collect(0:1:n_ph-1)
    # 2. Wiggly time (recomb):
    xdec = bg.x_grid[argmin(abs.( -ih.τ′ .* bg.ℋ .*bg.η .- 1))] #decoupling ish
    x_rc_f, n_rc = xdec, ie.Nᵧ₂ #100
    dx_rc = (x_rc_f-x_ph_f)/n_rc
    x_rc = x_ph_f .+ dx_rc * collect(1:1:n_rc)

    # 3. Post-recomb:
    n_pr = ie.Nᵧ₃ #50
    dx_pr = (bg.x_grid[end] -x_rc_f)/n_pr
    x_pr = x_rc_f .+ dx_pr* collect(1:1:n_pr )
    x_sparse = vcat(x_ph,x_rc,x_pr)
    return x_sparse
end


function boltsolve(ie::IE{T}, ode_alg=KenCarp4(); reltol=1e-6) where T #MD...
    x_grid = x_grid_ie(ie)
    xᵢ = first(x_grid)#ie.bg.x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=x_grid,#ie.bg.x_grid, 
                dense=false,
                )
    return sol
end


function boltsolve(ie::IEν{T}, ode_alg=KenCarp4(); reltol=1e-6) where T 
    xᵢ = first(ie.bg.x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=ie.bg.x_grid, dense=false,
                )
    return sol
end

struct ConformalIE{T<:Real,  H <: IE{T}, IT <: AbstractInterpolation{T}}
    ie::H
    η2x::IT
end

function boltsolve_conformal(confie::ConformalIE{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    ode_alg=KenCarp4(); reltol=1e-6) where T
    ie = confie.ie
    xᵢ = confie.η2x( ie.bg.η[1] ) #to be consistent
    u₀ = initial_conditions(xᵢ, ie)
    Mpcfac = bg.H₀*299792.458/100.
    prob = ODEProblem{true}(ie_conformal!, u₀, 
        (hierarchy.bg.η(hierarchy.bg.x_grid[1])*Mpcfac₀ , hierarchy.bg.η(hierarchy.bg.x_grid[end])*Mpcfac),
        confhierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,
    saveat=ie.bg.x_grid,
    dense=false
    )
    return sol
end

function ie_conformal!(du, u, confie::ConformalIE{T}, η) where T
    ie = confie.ie
    Mpcfac = ie.bg.H₀*299792.458/100.
    x = confie.η2x(η  / Mpcfac )
    ℋ = ie.bg.ℋ(x)
    ie!(du, u, ie, x)
    du .*= ℋ / Mpcfac  # account for dx/dη
    return nothing
end

function rsa_perts!(u, ie::IE{T},x) where T
    #redundant code for what we need to compute RSA perts in place in u
    k, ℓᵧ, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)
    ℓ_ν = ie.ℓ_ν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)

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


#---

function boltsolve_rsa(ie::IE{T}, ode_alg=KenCarp4(); reltol=1e-6) where T
    #call solve as usual first
    perturb = boltsolve(ie, reltol=reltol)
    x_grid = ie.bg.x_grid
    pertlen = 2(2+1)+(ie.ℓ_ν+1)+(ie.ℓ_mν+1)*ie.nq+5
    results=zeros(pertlen,length(x_grid))
    for i in 1:length(x_grid) results[:,i] = perturb(x_grid[i]) end
    #replace the late-time perts with RSA approx (assuming we don't change rsa switch)
    xrsa_hor = minimum(bg.x_grid[(@. k*bg.η .> 45)])
    xrsa_od = minimum(bg.x_grid[(@. -ih.τ′*bg.η*bg.ℋ .<5)])
    this_rsa_switch = max(xrsa_hor,xrsa_od)
    x_grid_rsa = x_grid[x_grid.>this_rsa_switch]
    results_rsa = results[:,x_grid.>this_rsa_switch]
    #(re)-compute the RSA perts so we can write them to the output vector
    for i in 1:length(x_grid_rsa) #inside here use regular unpack since single step
        rsa_perts!(view(results_rsa,:,i),ie,x_grid_rsa[i]) #to mutate need to use view...
    end
    results[:,x_grid.>this_rsa_switch] = results_rsa
    sol = results
    return sol
end

# basic Newtonian gauge: establish the order of perturbative variables in the ODE solve
function unpack(u, ie::IE{T, BasicNewtonian}) where T
    ℓ_ν =  ie.ℓ_ν
    ℓ_mν = ie.ℓ_mν #should be smaller than others
    nq = ie.nq
    ℓᵧ=2 
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    𝒩 = OffsetVector(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1)) , 0:ℓ_ν)  # indexed 0 through ℓ_ν
    ℳ = OffsetVector(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq )) , 0:(ℓ_mν+1)*nq -1)  # indexed 0 through ℓ_mν
    Φ, δ, v, δ_b, v_b = view(u, ((2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5)) #getting a little messy...
    return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
end

function unpack(u, ie::IEν{T, BasicNewtonian}) where T
    ℓᵧ =  ie.ℓ_γ
    ℓ_mν = ie.ℓ_mν #should be smaller than others
    nq = ie.nq
    ℓ_ν=2 
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    𝒩 = OffsetVector(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1)) , 0:ℓ_ν)  # indexed 0 through ℓ_ν
    ℳ = OffsetVector(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq )) , 0:(ℓ_mν+1)*nq -1)  # indexed 0 through ℓ_mν
    Φ, δ, v, δ_b, v_b = view(u, ((2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5)) #getting a little messy...
    return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
end

#FIXME this is probably terrible for performance
function ie_unpack(u, ie::IE{T, BasicNewtonian}) where T
    ℓ_ν =  ie.ℓ_ν
    ℓ_mν = ie.ℓ_mν #should be smaller than others
    nq = ie.nq
    Nᵧ = ie.Nᵧ₁+ie.Nᵧ₂+ie.Nᵧ₃ #ie.Nᵧ
    ℓᵧ=2
    #here u is the history of u over all ie timesteps (perlen,ie timesteps)
    #The perts below will be their histories over all ie timesteps as well
    #leading index will be pert index, trailing the time index
    Θ = OffsetArray(view(u, 1:(ℓᵧ+1),:), 0:ℓᵧ, 1:Nᵧ)  # indexed 0 through ℓᵧ, 1 through Nᵧ
    Θᵖ = OffsetArray(view(u, (ℓᵧ+2):(2ℓᵧ+2),:), 0:ℓᵧ, 1:Nᵧ)  # indexed 0 through ℓᵧ
    𝒩 = OffsetArray(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1),:) , 0:ℓ_ν, 1:Nᵧ)  # indexed 0 through ℓ_ν
    ℳ = OffsetArray(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq ),:) , 0:(ℓ_mν+1)*nq-1, 1:Nᵧ)  # indexed 0 through ℓ_mν
    # Φ, δ, v, δ_b, v_b = view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5, :) #getting a little messy...
	Φ, δ, v, δ_b, v_b = eachrow( view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5, :) ) #getting a little messy...
	return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b #perts over all ie timesteps
end

function ie_unpack(u, ie::IEν{T, BasicNewtonian}) where T
    ℓᵧ =  ie.ℓ_γ
    ℓ_mν = ie.ℓ_mν #should be smaller than others
    nq = ie.nq
    N_ν = ie.N_ν
    ℓ_ν=2
    Θ = OffsetArray(view(u, 1:(ℓᵧ+1),:), 0:ℓᵧ, 1:N_ν)  # indexed 0 through ℓᵧ, 1 through Nᵧ
    Θᵖ = OffsetArray(view(u, (ℓᵧ+2):(2ℓᵧ+2),:), 0:ℓᵧ, 1:N_ν)  # indexed 0 through ℓᵧ
    𝒩 = OffsetArray(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1),:) , 0:ℓ_ν, 1:N_ν)  # indexed 0 through ℓ_ν
    ℳ = OffsetArray(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq ),:) , 0:(ℓ_mν+1)*nq-1, 1:N_ν)  # indexed 0 through ℓ_mν
	Φ, δ, v, δ_b, v_b = eachrow( view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5, :) ) #getting a little messy...
	return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b #perts over all ie timesteps
end

function ie!(du, u, ie::IE{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih, ie.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)


    ℓ_ν = ie.ℓ_ν
    ℓ_mν =  ie.ℓ_mν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, ie)  # will be sweetened by .. syntax in 1.6
    Θ[2] = ie.sΘ2(x)# call the spline, update Θ₂ at top since we do not evolve it


    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
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

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)

    # neutrinos (massive, MB 57)
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
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
    rsa_on = false#(k*ηₓ > 45) &&  (-τₓ′*ηₓ*ℋₓ < 5)
    #*sqrt(H₀²)< 1) #is this ℋ or H0?
    if rsa_on
        # println("INSIDE RSA")
        #photons
        Θ[0] = Φ - ℋₓ/k *τₓ′ * v_b
        # Θ[1] = -2Φ′/k + (k^-2)*( τₓ′′ * v_b + τₓ′ * (ℋₓ*v_b - csb² *δ_b/k + k*Φ) )
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

        # manual zeroing to avoid saving garbage
        𝒩′[:] = zeros(ℓ_ν+1)
        Θ′[:] = zeros(ℓᵧ+1)
        Θᵖ′[:] = zeros(ℓᵧ+1)

    else
        #do usual ie
        # relativistic neutrinos (massless)
        𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
        𝒩′[1] = k/(3ℋₓ) * 𝒩[0] - 2*k/(3ℋₓ) *𝒩[2] + k/(3ℋₓ) *Ψ
        for ℓ in 2:(ℓ_ν-1)
            𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) * ( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1] )
        end
        #truncation (same between MB and Callin06/Dodelson)
        𝒩′[ℓ_ν] =  k / ℋₓ  * 𝒩[ℓ_ν-1] - (ℓ_ν+1)/(ℋₓ *ηₓ) *𝒩[ℓ_ν]


        # photons
        #Temp IE:
        # Θ[2] = IE_solve(∫Θ₂,xᵢ,x,Nᵧ) #how to choose xᵢ?

		#ℓ=0,1 DEs
        Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
        Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)

        # polarized photons
        #Polzn IE:
        # Π = IE_solve(∫Π,xᵢ,x,Nᵧ) #not doing the internal solve rn, try later
		Π = ie.sΠ(x) #call the spline
        Θᵖ[2] = Π - Θᵖ[0] - Θ[2]#get Θᵖ′[2] from Π again - this easy?
		#ℓ=0,1 DEs
        Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
        Θᵖ′[1] = k / (3ℋₓ) * Θᵖ[0] - 2k / (3ℋₓ) * Θᵖ[2] + τₓ′ * Θᵖ[1] #usual ie term but just for ℓ=1

    end
    #END RSA

    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end


function ie!(du, u, ie::IEν{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓ_ν, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih, ie.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)


    ℓᵧ = ie.ℓ_γ
    ℓ_mν =  ie.ℓ_mν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, ie)  # will be sweetened by .. syntax in 1.6
    # Θ[2] = ie.sΘ2(x)# call the spline, update Θ₂ at top since we do not evolve it
    𝒩[0] = ie.s𝒩₀(x)
    𝒩[2] = ie.s𝒩₂(x)


    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    # ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par)
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

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)

    # neutrinos (massive, MB 57)
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
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

    # relativistic neutrinos (massless)
    𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′
    𝒩′[1] = k/(3ℋₓ) * 𝒩[0] - 2*k/(3ℋₓ) *𝒩[2] + k/(3ℋₓ) *Ψ
    # for ℓ in 2:(ℓ_ν-1)
        # 𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) * ( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1] )
    # end
    #truncation (same between MB and Callin06/Dodelson)
    # 𝒩′[ℓ_ν] =  k / ℋₓ  * 𝒩[ℓ_ν-1] - (ℓ_ν+1)/(ℋₓ *ηₓ) *𝒩[ℓ_ν]


    # photons (hierarchy way)
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


    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end



#FIXME: don't need to copy all this code?
function initial_conditions(xᵢ, ie::IE{T, BasicNewtonian}) where T
    k, ℓᵧ, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih, ie.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    ℓ_ν = ie.ℓ_ν
    ℓ_mν =  ie.ℓ_mν
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ are mutable views (see unpack)

    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    #These get a 3/3 since massive neutrinos behave as massless at time of ICs
    Ω_ν =  7*(3/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
    f_ν = 1/(1 + 1/(7*(3/3)*par.N_ν/8 *(4/11)^(4/3)))

    # metric and matter perturbations
    Φ = 1.0
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) )

    #trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
   

    # ->apparently this does nothing TO BE CONSISTENT (even though this will give wrong ICs?)
    # Θ[2] = ie.sΘ2(xᵢ)# call the spline, update Θ₂ at top since we do not evolve it
    # Π = ie.sΠ(xᵢ) #call the spline
    
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    # TO BE CONSISTENT (even though this will give wrong ICs?)
    # Θᵖ[2] = Π - Θᵖ[0] - Θ[2]#get Θᵖ′[2] from Π again - this easy?

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
    for ℓ in 3:ℓ_ν
        𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #standard truncation
    end

    #massive neutrino hierarchy
    #It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
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

#FIXME this is a waste since the  only thing that changes is ℓ_ν vs ℓᵧ...
function initial_conditions(xᵢ, ie::IEν{T, BasicNewtonian}) where T
    k, ℓ_ν, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih, ie.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    ℓᵧ = ie.ℓ_γ
    ℓ_mν =  ie.ℓ_mν
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ are mutable views (see unpack)

    H₀²,aᵢ² = bg.H₀^2,exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    #These get a 3/3 since massive neutrinos behave as massless at time of ICs
    Ω_ν =  7*(3/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
    f_ν = 1/(1 + 1/(7*(3/3)*par.N_ν/8 *(4/11)^(4/3)))

    # metric and matter perturbations
    Φ = 1.0
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) )

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

    # neutrino hierarchy
    # we need xᵢ to be before neutrinos decouple, as always
    𝒩[0] = Θ[0]
    𝒩[1] = Θ[1]
    𝒩[2] = - (k^2 *ηₓ^2)/15 * 1 / (1 + 2/5 *f_ν) * Φ  / 2 #MB
    #FIXME^put the C here for consistency
    # for ℓ in 3:ℓ_ν
    #     𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #standard truncation
    # end

    #massive neutrino hierarchy
    #It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
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
#FIXME ignore source functions for now - nothing will need to change except struct arg

#---

# Spherical Bessel functions
j2(x) = (x > 0.01) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^3 : x^2 /15 - x^4 /210 + x^6 /7560
j2bx2(x)  =  (x > 0.06) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^5 : 1/15 - x^2 /210 + x^4 /7560 - x^6 /498960
j2′(x) = (x > 0.01) ? ( -x*(x^2 -9)*cos(x) + (4x^2 -9)*sin(x) ) / x^4 : 2x /15 - 2x^3 /105 + x^5 /1260
j2′′(x) = (x > 0.2) ? ( x*(5x^2 -36)*cos(x) + (x^4 - 17x^2 +36)*sin(x) ) / x^5 : 2/15 - 2x^2 /35 + x^4 /252 - x^6 /8910
R2(x) = (x > 0.2) ? -( j2(x) + 3j2′′(x) ) / 2 : -1/5 + 11x^2 /210 -x^4 /280 +17x^4 /166320

#IE helper function
function get_perts(u,ie::IE{T},x) where T
    k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᴾ, 𝒩, ℳ, Φ, δ, _, δ_b, v_b = unpack(u, ie)  # the Θ, Θᴾ, 𝒩 are views (see unpack)

    ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
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

    #update with splines
    Θ[2] = ie.sΘ2(x)
    Π = ie.sΠ(x)
    Θᴾ[2] = Π - Θᴾ[0] - Θ[2]
    return Φ′,Ψ,Θ[0],Π,v_b
end

#Kamionkowski weights
function Ws(xⱼ,xⱼ₊₁,τ,xᵢ)
    ϵτ = 1e-3 #if g is changing rapidly (τ′*dx>ϵτ), use g-aware weights #FIXME is the switch really necessary?
    dτ = -(τ(xⱼ₊₁) - τ(xⱼ))  #NB this is backwards from Kamionkowski since he does ``from 1''
    expτ = exp(-dτ)
    expj = exp( -( -τ(xᵢ) + τ(xⱼ₊₁) ) ) #NB ditto above
    τfac = (1 - (1+dτ)*expτ)/ dτ
    Wⱼ⁺ = expj* (  (dτ >ϵτ)  ? 1-expτ - τfac : dτ/2  )
    Wⱼ =  expj* (  (dτ >ϵτ)  ? τfac : dτ/2  )
    return Wⱼ,Wⱼ⁺
end

# KERNELS
function _IΘ2(x, x′,k,
    Π, Θ0, v_b, Φ′, Ψ,
    ih, bg) #for testing
    τ′,η = ih.τ′,bg.η #all splines of x
    y = k*( η(x)-η(x′) )#Bessel argument
    IΘ2 = ( Θ0 - Φ′/ (-τ′(x′))  )*j2(y) - ( v_b   - ( k/bg.ℋ(x′) )*Ψ / (-τ′(x′)) )*j2′(y)  - Π*R2(y) / 2 
    return IΘ2
end

function _IΠ(x, x′,k, Π, bg)
    η = bg.η #all splines of x
    y = k*( η(x)-η(x′) )#Bessel argument
    IE2 = j2bx2(y)*Π
    IΠ = 9IE2
    return IΠ
end
function g_weight_trapz_ie(i,x_grid,ie::IE{T},Φ′,Ψ,Θ₀,Π,v_b) where T
    τ = ie.ih.τ
    k = ie.k
    xᵢ = x_grid[i]
    Θ2ᵢ, Πᵢ = 0,0 
    Wⱼ,Wⱼ⁺ = 0, 0 
    for j in 1:i-2
        xⱼ,xⱼ₊₁ = x_grid[j], x_grid[j+1]
        Wⱼ,Wⱼ⁺ = Ws(xⱼ,xⱼ₊₁,τ,xᵢ) #passing xᵢ for now but could update later externally...
        #TODO if we want to compute weights once for all i,j and save them we can?
        #Implicit weighting scheme at each timestep
        Θ2ᵢ += (_IΘ2(xᵢ,xⱼ₊₁,k,Π[j+1],Θ₀[j+1],v_b[j+1],Φ′[j+1],Ψ[j+1],ie.ih,ie.bg)*Wⱼ⁺
               + _IΘ2(xᵢ,xⱼ,k,Π[j],Θ₀[j],v_b[j],Φ′[j],Ψ[j],ie.ih,ie.bg)*Wⱼ)
        Πᵢ += ( _IΠ(xᵢ,xⱼ₊₁,k,Π[j+1],ie.bg)*Wⱼ⁺
               + _IΠ(xᵢ,xⱼ,k,Π[j],ie.bg)*Wⱼ)
    end
    #Handle final sub-timestep j = i-1 (pull out final loop iteration)
    xⱼ,xⱼ₊₁ = x_grid[i-1], xᵢ
    Wⱼ,Wⱼ⁺ = Ws(xⱼ,xⱼ₊₁,τ,xᵢ) #passing xᵢ for now but could update later externally...
    Θ2ᵢ += (_IΘ2(xᵢ,x_grid[i],k,0.,Θ₀[i],v_b[i],Φ′[i],Ψ[i],ie.ih,ie.bg)*Wⱼ⁺
           + _IΘ2(xᵢ,x_grid[i-1],k,Π[i-1],Θ₀[i-1],v_b[i-1],Φ′[i-1],Ψ[i-1],ie.ih,ie.bg)*Wⱼ)
    Πᵢ += _IΠ(xᵢ,x_grid[i-1],k,Π[i-1],ie.bg)*Wⱼ
    #Kamionkowski integration scheme for handling x′ = x at each x (the implicit timestep)
    Πᵢ = (Πᵢ + Θ2ᵢ) / ( 1 - 7/10 * Wⱼ⁺)
    Θ2ᵢ = Θ2ᵢ + ( i <length(x_grid) ? Π[i+1] : 0. )/10 * Wⱼ⁺ #if i+1>length(x_grid), return 0 for oob array
    return Θ2ᵢ,Πᵢ
end

#FIXME? consolidate on interpolator -> pass an interpolator rather than the ingredients?


function iterate(Θ₂_km1,Π_km1, 𝕡::CosmoParams{T}, bg, ih, k, 
    Nᵧ₁,Nᵧ₂,Nᵧ₃,xgi,
    ℓ_ν, ℓ_mν, n_q,reltol) where T
    Θ₂_k,Π_k = zero(Θ₂_km1),zero(Π_km1) #maybe pre-allocate these (and below)
    ie_k = IE(BasicNewtonian(), 𝕡, bg, ih, k,
            linear_interpolation(xgi,Θ₂_km1),
            linear_interpolation(xgi,Π_km1),
            Nᵧ₁,Nᵧ₂,Nᵧ₃,
            ℓ_ν, ℓ_mν, n_q)
    u_all_k = boltsolve(ie_k; reltol=reltol)
    N = length(xgi)
    Φ′,Ψ,Θ₀,Π,v_b = zeros(N),zeros(N),zeros(N),zeros(N),zeros(N)
    for (j,u) in enumerate( eachcol(u_all_k) )
            Φ′[j],Ψ[j],Θ₀[j],Π[j],v_b[j] = get_perts(u,ie_k,xgi[j])
    end
    for i in 3:length(xgi)
            Θ₂_k[i],Π_k[i] = g_weight_trapz_ie(i,xgi,ie_k,Φ′,Ψ,Θ₀,Π,v_b)
    end
    return Θ₂_k,Π_k,u_all_k
end