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
    Nᵧ::Int #can't we just use existing x grid for this?
    ℓ_ν::Int
    ℓ_mν::Int
    nq::Int
end

IE(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real,
    sΘ2::AbstractInterpolation,sΠ::AbstractInterpolation,
	#^FIXME: Are these right?? I dropped the {T,1} since T is not known here to get it to compile
    Nᵧ=400, ℓ_ν=8, ℓ_mν=10, nq=15
    ) = IE(integrator, par, bg, ih, k, sΘ2, sΠ, Nᵧ, ℓ_ν,ℓ_mν, nq)


function iesolve(ie,u,Θ₂,Π)
	#ie is the existing integrator
	#u is history of perts
	xx = ie.bg.x_grid
	# Θ₂,Π = zeros(ie.Nᵧ), zeros(ie.Nᵧ) #local arrays #TODO maybe make ie attributes?
    for i in 2:length(xx)
		Θ₂[i],Π[i] = g_weight_trapz_ie(xx[i],ie,u[i])
	end
end

function itersolve(ie)
    # initialize ansatz - start with zero
	Θ₂,Π = zeros(ie.Nᵧ),zeros(ie.Nᵧ)
    # start picard iteration

        # update splines
		ie.Θ₂ = spline(Θ₂, ie.bg.x_grid)
		ie.Π =  spline(Π,  ie.bg.x_grid)
        # solve odes
		perturb = boltsolve_rsa(ie)

        # solve ie to get Θ₂, Π - Picard step with weights for coupling is here
		iesolve(ie,perturb.u,Θ₂,Π)

end


function boltsolve(ie::IE{T}, ode_alg=KenCarp4(); reltol=1e-6) where T #MD...
    xᵢ = first(ie.bg.x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=ie.bg.x_grid, dense=false,
                )
    return sol
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

function boltsolve_rsa(ie::IE{T}, ode_alg=KenCarp4(); reltol=1e-6) where T
    #call solve as usual first
    perturb = boltsolve(ie, reltol=reltol)
    x_grid = ie.bg.x_grid
    pertlen = 2(2+1)+(ie.ℓ_ν+1)+(ie.ℓ_mν+1)*ie.nq+5
    results=zeros(pertlen,length(x_grid))
    for i in 1:length(x_grid) results[:,i] = perturb(x_grid[i]) end
    #replace the late-time perts with RSA approx (assuming we don't change rsa switch)
    this_rsa_switch = x_grid[argmin(abs.(ie.k .* ie.bg.η.(x_grid) .- 45))]
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

#FIXME this is probably terrible for performance
function ie_unpack(u, ie::IE{T, BasicNewtonian}) where T
    ℓ_ν =  ie.ℓ_ν
    ℓ_mν = ie.ℓ_mν #should be smaller than others
    nq = ie.nq
    Nᵧ = ie.Nᵧ
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
    # println("k condition ", k*ηₓ)
    # println("tau condition ", -5τₓ′*ηₓ*ℋₓ)
    # if (k*ηₓ > 45) println("k condition satisfied") end
    # if -5τₓ′*ηₓ*sqrt(H₀²)< 1 println("tau condition satisfied") end
    rsa_on = false#(k*ηₓ > 45) &&  (-5τₓ′*ηₓ*ℋₓ<1)
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

#FIXME need to import bessel functions somewhere?

# The RHSs of the IEs
# function ∫Θ₂(,,)
#
# end
#
# function ∫Π(,,)
#
# end
#
# # Volterra solver
# function IE_solve(∫f,N)
#
# end



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

#FIXME ignore source functions for now - nothing will need to change except struct arg
