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
    # N_ν::Int we don't need this I think
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
    # N_ν=400, 
    ℓ_γ=8, ℓ_mν=10, nq=15
    ) = IEν(integrator, par, bg, ih, k, s𝒩₀, s𝒩₂, N_ν, ℓ_γ,ℓ_mν, nq)


struct ConformalIE{T<:Real,  H <: IE{T}, IT <: AbstractInterpolation{T}}
        ie::H
        η2x::IT
    end

#lazy copy for now...
struct ConformalIEν{T<:Real,  H <: IEν{T}, IT <: AbstractInterpolation{T}}
        ie::H
        η2x::IT
    end
    


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
    xhor = bg.x_grid[argmin(abs.(k .* bg.η .- 2π))] #horizon crossing ish
    x_ph_i, x_ph_f, n_ph = bg.x_grid[1], xhor, ie.Nᵧ₁ #10.
    dx_ph = (x_ph_f-x_ph_i)/(n_ph-1)
    # x_ph = -20.:dx_ph:x_ph_f
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

function η_grid_ie(ie,η2x,N) 
    dx_η = (ie.bg.η[end] - ie.bg.η[1])/(N-1)
    #FIXME for loop version doesn't work...must be off by 1 or something
    # xx = zeros(N)
    # for i in 0:N-1
    #     xx[i+1] = η2x(ie.bg.η[1] + dx_η*i)
    # end
    ηs = ie.bg.η[1] .+ dx_η* collect(0:1:N-1 )
    # return xx
    return η2x.(ηs)
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
                saveat=ie.bg.x_grid, dense=false, #FIXME
                )
    return sol
end


function boltsolve_conformal(confie::ConformalIE{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    x_grid = η_grid_ie(ie,η2x,2048) #this is overkill/unoptomized but just to have something that decently agrees...
    xᵢ = first(x_grid) #to be consistent
    # xᵢ = confie.η2x( ie.bg.η[1] ) 
    u₀ = initial_conditions(xᵢ, ie)
    Mpcfac = ie.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(ie_conformal!, u₀, 
        (ie.bg.η(xᵢ)*Mpcfac, ie.bg.η(x_grid[end])*Mpcfac),
        confie)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=ie.bg.η(x_grid)*Mpcfac,
                dense=false
                )
    return sol
end

function boltsolve_conformal(confie::ConformalIEν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    # x_grid = η_grid_ie(ie,η2x,2048) #why are we doing this? do something instead similar to perturbations.jl
    # xᵢ = first(x_grid) #to be consistent
    xᵢ = η2x( ie.bg.η[1] ) 
    Mpcfac = ie.bg.H₀*299792.458/100.
    # xᵢ = η2x( 1.0/Mpcfac ) 
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie_conformal!, u₀, 
        (ie.bg.η[1]*Mpcfac, ie.bg.η[end]*Mpcfac),
        confie)
    sol = solve(prob, ode_alg, reltol=reltol,
                # saveat=ie.bg.η(x_grid)*Mpcfac,
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

function ie_conformal!(du, u, confie::ConformalIEν{T}, η) where T
    ie = confie.ie
    Mpcfac = ie.bg.H₀*299792.458/100.
    x = confie.η2x(η  / Mpcfac )
    ℋ = ie.bg.ℋ(x)
    ie!(du, u, ie, x)
    du .*= ℋ / Mpcfac  # account for dx/dη
    return nothing
end

function itersolve_conformal(Nₖ::Int,confie::ConformalIE{T};reltol=1e-6) where T
    ie_0, η2x = confie.ie, confie.η2x
    #All we have to do is change the time points to be equispaced in \eta
    x_grid = η_grid_ie(ie_0,η2x,2048) #1000 is not great, 2048 is good not perfect, leave for now
    # println("xgrids: ",x_grid[1],", ",x_grid[end])
    Θ₂,Π =  zeros(length(x_grid)),zeros(length(x_grid)) #initialize to zero (for now)
    pertlen = 2(2+1)+(ie_0.ℓ_ν+1)+(ie_0.ℓ_mν+1)*ie_0.nq+5
    u_all = zeros(pertlen,length(x_grid))
    #FIXME - make In-place?
    for k in 1:Nₖ 
            Θ₂,Π,u_all = iterate_conformal(Θ₂,Π, ie_0.par, ie_0.bg, ie_0.ih, ie_0.k,  
                                    ie_0.Nᵧ₁,ie_0.Nᵧ₂,ie_0.Nᵧ₃, 
                                    x_grid, ie_0.ℓ_ν, ie_0.ℓ_mν, ie_0.nq, 
                                    reltol,
                                    η2x)
    end
    return u_all
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

# Dead code
# function ie_unpack(u, ie::IEν{T, BasicNewtonian}) where T
#     ℓᵧ =  ie.ℓ_γ
#     ℓ_mν = ie.ℓ_mν #should be smaller than others
#     nq = ie.nq
#     N_ν = ie.N_ν
#     ℓ_ν=2
#     Θ = OffsetArray(view(u, 1:(ℓᵧ+1),:), 0:ℓᵧ, 1:N_ν)  # indexed 0 through ℓᵧ, 1 through Nᵧ
#     Θᵖ = OffsetArray(view(u, (ℓᵧ+2):(2ℓᵧ+2),:), 0:ℓᵧ, 1:N_ν)  # indexed 0 through ℓᵧ
#     𝒩 = OffsetArray(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1),:) , 0:ℓ_ν, 1:N_ν)  # indexed 0 through ℓ_ν
#     ℳ = OffsetArray(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq ),:) , 0:(ℓ_mν+1)*nq-1, 1:N_ν)  # indexed 0 through ℓ_mν
# 	Φ, δ, v, δ_b, v_b = eachrow( view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5, :) ) #getting a little messy...
# 	return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b #perts over all ie timesteps
# end

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
    
    Mpcfac = ie.bg.H₀*299792.458/100.
    # if ηₓ*Mpcfac >= 1.0 #overwrite the neutrino perts if sufficiently late
    𝒩[0] = ie.s𝒩₀(x) #FIXME this sucks, need a ctime ie!
    𝒩[2] = ie.s𝒩₂(x)
    # end

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
    # println("type Phi: ", typeof(Φ′))
    # println("type Phi2 : ", typeof(Φ′) <: ForwardDiff.Dual)

    # if ((x<=-19.99 || x>=-0.01) &&  ~(typeof(Φ′) <: ForwardDiff.Dual))
    #     println("x = ", x)
    #     println("Φ′ = ", Φ′)
    #     println("𝒩[0] after = ", 𝒩[0])
    #     println("Θ[0] = ", Θ[0])
    #     println("𝒩[2] after  = ", 𝒩[2])
    #     println("Ψ = ", Ψ)
    #     println("Ψ components: Θ₂ = $(Θ[2]), 𝒩₂ = $(𝒩[2]), σℳ = $(σℳ)")
    # end

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
    #use truncation expression since we don't evolve octopole
    # if ηₓ*Mpcfac < 1.0  #if early, need to actually evolve quadrupole
    #     𝒩′[2] =  k / ℋₓ  * 𝒩[1] - 3/(ℋₓ *ηₓ) *𝒩[2]
    # end #if later, just don't evolve this, result will be junk
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

    #update with splines #FIXME does this actually do anything?
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

#FIXME? consolidate on interpolator -> pass an interpolator rather than the pert ingredients?


function iterate(Θ₂_km1,Π_km1, 𝕡::CosmoParams{T}, bg, ih, k, 
    Nᵧ₁,Nᵧ₂,Nᵧ₃,xgi,
    ℓ_ν, ℓ_mν, n_q,reltol) where T
    Θ₂_k,Π_k = zero(Θ₂_km1),zero(Π_km1) #FIXME pre-allocate these (and below)
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

function iterate_conformal(Θ₂_km1,Π_km1, 𝕡::CosmoParams{T}, bg, ih, k, 
    Nᵧ₁,Nᵧ₂,Nᵧ₃,xgi,
    ℓ_ν, ℓ_mν, n_q,reltol,
    η2x) where T
    Θ₂_k,Π_k = zero(Θ₂_km1),zero(Π_km1) #FIXME pre-allocate these (and below)
    ie_k = IE(BasicNewtonian(), 𝕡, bg, ih, k,
            linear_interpolation(xgi,Θ₂_km1),
            linear_interpolation(xgi,Π_km1),
            Nᵧ₁,Nᵧ₂,Nᵧ₃,
            ℓ_ν, ℓ_mν, n_q)
    ie_k_conf = ConformalIE(ie_k,η2x);
    u_all_k = boltsolve_conformal(ie_k_conf; reltol=reltol)
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

# ------------------------------
# FFT Iteration functions

# # /// IC Free streaming ///
# # Relevant Bessel functions (ℓ=0,1,2)
# #ℓ=0
# j0(x) = (x > 0.01) ? sin(x)/x : 1 - x^2 /6 + x^4 /120 - x^6 /5040
# j0′(x) = -j1(x)
# #ℓ=1
# j1(x) =  (x > 0.01) ?  (sin(x) - x*cos(x))/x^2 : x/3 - x^3 /30 + x^5 /840
# R1(x) =  (x > 0.01) ? j1(x) - 3j2(x)/x : 2x/15 - 2x^3 /105 + x^5 /1260
# #ℓ=2
# j2(x) = (x > 0.01) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^3 : x^2 /15 - x^4 /210 + x^6 /7560
# j2′(x) = (x > 0.01) ? ( -x*(x^2 -9)*cos(x) + (4x^2 -9)*sin(x) ) / x^4 : 2x /15 - 2x^3 /105 + x^5 /1260
# j2′′(x) = (x > 0.2) ? ( x*(5x^2 -36)*cos(x) + (x^4 - 17x^2 +36)*sin(x) ) / x^5 : 2/15 - 2x^2 /35 + x^4 /252 - x^6 /8910
# R2(x) = (x > 0.2) ? -( j2(x) + 3j2′′(x) ) / 2 : -1/5 + 11x^2 /210 -x^4 /280 +17x^4 /166320
# # The W coupling kernel (sum truncated at ℓ=2)
# W00(x) = j0(x)
# W01(x) = j1(x)
# W02(x) = j2(x)
# W21(x) = -R1(x)
# W22(x) = -R2(x)
# function Wsum(x,𝒳ᵢ₀,𝒳ᵢ₁,𝒳ᵢ₂)
#     𝒳ₛ₀ = W00(x)*𝒳ᵢ₀ - 3W01(x)*𝒳ᵢ₁ + 5W02(x)*𝒳ᵢ₂  #ℓ=0 ( use the subscript ₛ for streaming, this is the "free-streaming" piece)
#     𝒳ₛ₂ = W02(x)*𝒳ᵢ₀ - 3W21(x)*𝒳ᵢ₁ + 5W22(x)*𝒳ᵢ₂ #ℓ=2
#     return 𝒳ₛ₀, 𝒳ₛ₂
# end

# function get_Φ′_Ψ(u,hierarchy::Hierarchy{T},x) where T
#     #TODO: can streamline hierarchy and source funcs with this helper function also
#     k, par, bg, nq = hierarchy.k, hierarchy.par, hierarchy.bg,hierarchy.nq
#     Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
#     ℋₓ =  bg.ℋ(x)
#     a = x2a(x)
#     Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
#     Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
#     ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
#     Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
#                                   Ω_ν * 𝒩[2]
#                                   + σℳ / bg.ρ_crit /4
#                                   )
#     Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
#         Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
#         + 4Ω_r * a^(-2) * Θ[0]
#         + 4Ω_ν * a^(-2) * 𝒩[0]
#         + a^(-2) * ρℳ / bg.ρ_crit
#         )
#     return Φ′,Ψ
# end

# # Get the Φ' and Ψ (copy function in ie file) from hierarchy
# function get_Φ′_Ψ(u,ie::IEν{T},x) where T
#     #TODO: can streamline hierarchy and source funcs with this helper function also
#     k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
#     Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
#     ℋₓ =  bg.ℋ(x)
#     a = x2a(x)
#     Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
#     Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
#     𝒩[0] = ie.s𝒩₀(x)
#     𝒩[2] = ie.s𝒩₂(x)#WHY DO WE NEED THIS HERE BUT NOT IN PHOTONS? AND NOT FOR MONO?
#     ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
#     Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
#                                   Ω_ν * 𝒩[2]
#                                   + σℳ / bg.ρ_crit /4
#                                   )
#     Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
#         Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
#         + 4Ω_r * a^(-2) * Θ[0]
#         + 4Ω_ν * a^(-2) * 𝒩[0]
#         + a^(-2) * ρℳ / bg.ρ_crit
#         )
#     return Φ′,Ψ
# end

# function fft_funcs(x, y, Φ′,Ψ, k,ℋ,q,m,𝕡)
#     ϵ = (q^2 .+ exp.(2x)*m^2 ).^(1/2) #sqrt syntax doesn't work w/ bcast but put it back when undo bcast...
#     q̃ = ϵ/q #convenience notation
#     G₀ = ℋ .* q̃/k .* Φ′ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #for integrating in y #
#     G₁ = -q̃.^2 .* Ψ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #
#     K₀₀ = j0.(y) #1st index is ℓ 2nd index is derivative order
#     K₀₁ = j0′.(y)
#     K₂₀ = j2.(y) #
#     K₂₁ = j2′.(y) #
#     return G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁
# end

# function fft_integral(x, y,Φ′,Ψ,k,ℋ,q,m,𝕡,M) # for massive or massless neutrinos (𝒳=𝒩,ℳ)
#     dy = y[2]-y[1]
#     #  all ffts are performed in this function
#     G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁ = fft_funcs(x,y, Φ′,Ψ, k,ℋ,q,m,𝕡) #
#     # zero-pad the signals so convolution is not circular
#     G₀,G₁ = [G₀; zeros(M-1)],[G₁; zeros(M-1)]
#     K₀₀,K₀₁,K₂₀,K₂₁ = [K₀₀; zeros(M-1)],[K₀₁; zeros(M-1)],[K₂₀; zeros(M-1)],[K₂₁; zeros(M-1)] #
#     # FFT the Gs, Ks
#     G̃₀,G̃₁ = fft(G₀),fft(G₁)
#     K̃₀₀, K̃₀₁, K̃₂₀, K̃₂₁ = fft(K₀₀),fft(K₀₁),fft(K₂₀),fft(K₂₁)#
#     # Convolution theorem (iFFT pointwise product)
#     𝒳₀ₓ = ifft(G̃₀.*K̃₀₀ .+ G̃₁.*K̃₀₁)[1:M]*dy 
#     𝒳₂ₓ = ifft(G̃₀.*K̃₂₀ .+ G̃₁.*K̃₂₁)[1:M]*dy 
#     return 𝒳₀ₓ,𝒳₂ₓ
# end

# function fft_ie(ie,perturb,M,m,q,i_q,u₀,x_grid)
#     𝕡,bg,k,nq = ie_0.par,ie_0.bg,ie_0.k,ie.nq
#     # Set up the "neutrino horizon" and FFT abscissas
#     χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
#     yyx = k.* (χνs .- χνs[1])
#     dy=(yyx[end]-yyx[1])/(M-1)
#     yy = yyx[1]:dy:yyx[end]
#     invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
#     # Get metric sources
#     Φ′,Ψ = zeros(M),zeros(M)
#     for j in 1:M
#         Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
#     end
#     _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
#     if m==0 
#         𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
#     else
#         𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
#     end 
#     # Compute the new perts via FFT
#     𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
#     # Put it all together
#     𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
#     𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
#     return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)
# end

# function fft_ie_c(ie,perturb,M,m,q,i_q,u₀,x_grid)
#     𝕡,bg,k,nq = ie_0.par,ie_0.bg,ie_0.k,ie.nq
#     # Set up the "neutrino horizon" and FFT abscissas
#     χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
#     yyx = k.* (χνs .- χνs[1])
#     dy=(yyx[end]-yyx[1])/(M-1)
#     yy = yyx[1]:dy:yyx[end]
#     invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
#     # Get metric sources
#     Φ′,Ψ = zeros(M),zeros(M)
#     for j in 1:M
#         Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
#     end
#     _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
#     if m==0 
#         𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
#     else
#         𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
#     end 
#     # Compute the new perts via FFT
#     𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
#     # Put it all together
#     𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
#     𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
#     return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)#,
# end

# function h_boltsolve_flex(hierarchy::Hierarchy{T},  x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T
#     prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (x_ini , x_fin), hierarchy)
#     sol = solve(prob, ode_alg, reltol=reltol,
#                 dense=false,
#                 )
#     return sol
# end

# function boltsolve_flex(ie::IEν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
#     prob = ODEProblem{true}(Bolt.ie!, u₀, (x_ini , x_fin), ie)
#     sol = solve(prob, ode_alg, reltol=reltol,
#                 saveat=ie.bg.x_grid, dense=false, #FIXME
#                 )
#     return sol
# end

# function h_boltsolve_conformal_flex(confhierarchy::ConformalHierarchy{T},#FIXME we do't need this? {Hierarchy{T},AbstractInterpolation{T}},
#     η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
#     hierarchy = confhierarchy.hierarchy
#     Mpcfac = hierarchy.bg.H₀*299792.458/100.
#     prob = ODEProblem{true}(Bolt.hierarchy_conformal!, u₀, 
#                             (η_ini*Mpcfac , η_fin*Mpcfac),
#                             confhierarchy)
#     sol = solve(prob, ode_alg, reltol=reltol,
#     dense=false
#     )
#     return sol
# end

# function boltsolve_conformal_flex(confie::ConformalIEν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
#     η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
#     ie,η2x = confie.ie,confie.η2x
#     Mpcfac = ie.bg.H₀*299792.458/100.
#     prob = ODEProblem{true}(Bolt.ie_conformal!, u₀, 
#                             (η_ini*Mpcfac, η_fin*Mpcfac),
#                             confie)
#     sol = solve(prob, ode_alg, reltol=reltol,
#     dense=false
#     )
#     return sol
# end


# function iterate_fft(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
#     M, reltol,x_ini, x_fin,u0) where T
#     𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
#     ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
#                     𝒩₀_km1, 𝒩₂_km1,
#                     ℓᵧ, ℓ_mν, n_q)
#     perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
#     xx,𝒩₀_k,𝒩₂_k = fft_ie(ie_k_late,perturb_k_late,M,0.,1.,0,
#                         u0,perturb_k_late.t) #This is for massless only 
#     return xx,𝒩₀_k,𝒩₂_k,perturb_k_late
# end

# function iterate_fft_c(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
#     M, reltol,η_ini, η_fin,u0) where T
#     𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
#     ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
#                     𝒩₀_km1, 𝒩₂_km1,
#                     ℓᵧ, ℓ_mν, n_q)
#     ie_k_conf_late_c = ConformalIEν(ie_k_late,η2x);
#     perturb_k_late_c = boltsolve_conformal_flex(ie_k_conf_late_c, η_ini, η_fin, u0; reltol=reltol)
#     xx,𝒩₀_k,𝒩₂_k = fft_ie_c(ie_k_conf_late_c.ie,perturb_k_late_c,M,0.,1.,0,
#                         u0,η2x(perturb_k_late_c.t/Mpcfac)) #This is for massless only 
#     return xx,𝒩₀_k,𝒩₂_k,perturb_k_late_c
# end

# #---------------------------------#
# # Itersolves
# #---------------------------------#
# function itersolve_fft(Nₖ::Int,ie_0::IEν{T},M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
#     𝒩₀_0,𝒩₂_0 = ie_0.s𝒩₀,ie_0.s𝒩₂
#     𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
#     perturb_k = nothing
#     xx_k = nothing
#     for k in 1:Nₖ
#         xx_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft(𝒩₀_k,𝒩₂_k,ie_0.par,ie_0.bg,ie_0.ih,
#                                    ie_0.k,ie_0.ℓ_γ,ie_0.ℓ_mν,ie_0.nq,
#                                    M,reltol,x_ini,x_fin,u0)
#     end
#     return xx_k, 𝒩₀_k,𝒩₂_k,perturb_k
# end
# #ctime version
# function itersolve_fft(Nₖ::Int,ie_0_c::ConformalIEν{T},M::Int,η_ini, η_fin,u0;reltol=1e-6) where T
#     𝒩₀_0,𝒩₂_0 = ie_0_c.ie.s𝒩₀,ie_0_c.ie.s𝒩₂
#     𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
#     perturb_k = nothing
#     ηη_k = nothing
#     for k in 1:Nₖ
#         ηη_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft_c(𝒩₀_k,𝒩₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
#                                                ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.ℓ_mν,ie_0_c.ie.nq,M,reltol,
#                                                η_ini, η_fin,u0)
#     end
#     return ηη_k,𝒩₀_k,𝒩₂_k,perturb_k
# end

# # Helper functon for switch
# function get_switch_u0(η,hierarchy_conf) #Input is η of the switch
#     # switch_idx=740 #<- the switch idx for η=1.0ish
#     hierarchy,bg = hierarchy_conf.hierarchy,hierarchy_conf.hierarchy.bg
#     Mpcfac = bg.H₀*299792.458/100.
#     switch_idx = argmin(abs.(bg.η*Mpcfac .-η)) #for now we use the bg to find the switch
#     #solve the split ode
#     ℓᵧ,ℓ_ν,n_q = hierarchy.ℓᵧ,hierarchy.ℓ_ν, hierarchy.nq
#     pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
#     # \/ we want to report this timing to get a full picture of total time (early+late)
#     sol_early_c = h_boltsolve_conformal_flex(hierarchy_conf, bg.η[1], bg.η[switch_idx],  initial_conditions(bg.x_grid[1], hierarchy));
#     # Get the new initial conditions
#     u0_ie_c = zeros(2(ℓᵧ+1) + (2+1) + (ℓ_mν+1)*n_q + 5);
#     for i in  1:2(ℓᵧ+1)+(2+1)
#         u0_ie_c[i] = sol_early_c.u[end][i]
#     end
#     for i in  2(ℓᵧ+1)+(ℓ_ν+1)+1:pertlen
#         down_shift = i-(ℓ_ν-2)
#         u0_ie_c[down_shift] = sol_early_c.u[end][i]
#     end
#     return u0_ie_c
# end


# ------------------------------
# Unused
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
