# these types and functions integrate the Boltzmann hierarchy through time


#jms FIXME: In the great cleanup, I will:
# 1. merge all the different IE structs together into one big one
# 2. drop all the unnecessary functions that integrate in x rather than conformal time
# 3. turn to the issue of making ctime integration not some hacky thing

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
    ℓ_γ::Int
    ℓ_mν::Int
    nq::Int
end

#TODO the issue here is that you can't call T, IT{T}, AIT{IT}
# i.e. we can go one level but not 2, at least not without diff syntax
#trying a new example below...
# based  onn 
# https://stackoverflow.com/questions/25490364/method-will-not-match-with-nested-type-restrictions
# struct AI1{T <:Real, V <: AbstractInterpolation{T,1}}
#     AT::AbstractArray{V,1}
# end
# abstract type AbstractIE{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
#                          BG<:AbstractBackground, IH<:AbstractIonizationHistory, 
#                          Tk<:Real, #AT<:AbstractArray{T,1},
#                         IT<:AbstractInterpolation{T,1}, 
#                         } end
    # struct IEallν{T,PI,CP,BG,IH,Tk,IT,AIT
    #     } <: AbstractIE{T,PI,CP,BG,IH,Tk,IT,AIT
    #                      }
struct IEallν{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
        BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real,# AT<:AbstractArray{T,1},
        IT<:AbstractInterpolation{T,1}
    }
    integrator::PI
    par::CP
    bg::BG
    ih::IH
    k::Tk
    # sx::Array{T,1}
    # a𝒳₀:: AbstractArray{AT,1}
    # a𝒳₂::AbstractArray{AT,1}
    s𝒳₀::AbstractArray{IT,1}
    s𝒳₂::AbstractArray{IT,1}
    ℓ_γ::Int
    nq::Int
end

struct IEγν{T<:Real, PI<:PerturbationIntegrator, CP<:AbstractCosmoParams{T},
    BG<:AbstractBackground, IH<:AbstractIonizationHistory, Tk<:Real,
    IT<:AbstractInterpolation{T,1}
    }
    integrator::PI
    par::CP
    bg::BG
    ih::IH
    k::Tk
    sΘ2::IT #This is kept separate from the neutrino interpolators for convenience rn, but need not be
    sΠ::IT
    s𝒳₀::AbstractArray{IT,1}
    s𝒳₂::AbstractArray{IT,1}
    Nᵧ₁::Int #pre-entry
    Nᵧ₂::Int #recomb
    Nᵧ₃::Int #post-recomb
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
    ℓ_γ=8, ℓ_mν=10, nq=15
    ) = IEν(integrator, par, bg, ih, k, s𝒩₀, s𝒩₂, ℓ_γ,ℓ_mν, nq)

IEallν(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real,
    s𝒳₀::AbstractArray,s𝒳₂::AbstractArray,
    ℓ_γ=8, nq=15
    ) = IEallν(integrator, par, bg, ih, k, 
                s𝒳₀, s𝒳₂, 
                ℓ_γ, nq)

IEγν(integrator::PerturbationIntegrator, par::AbstractCosmoParams, bg::AbstractBackground,
    ih::AbstractIonizationHistory, k::Real,
    sΘ2::AbstractInterpolation,sΠ::AbstractInterpolation,
    s𝒳₀::AbstractArray,s𝒳₂::AbstractArray,
    Nᵧ₁=10, Nᵧ₂=100, Nᵧ₃=50, nq=15
    ) = IEγν(integrator, par, bg, ih, k, 
            sΘ2, sΠ,
            s𝒳₀, s𝒳₂, 
            Nᵧ₁,Nᵧ₂, Nᵧ₃, nq)

struct ConformalIE{T<:Real,  H <: IE{T}, IT <: AbstractInterpolation{T}}
        ie::H
        η2x::IT
    end

#lazy copy for now...
struct ConformalIEν{T<:Real,  H <: IEν{T}, IT <: AbstractInterpolation{T}}
        ie::H
        η2x::IT
    end
    
struct ConformalIEallν{T<:Real,  H <: IEallν{T}, IT <: AbstractInterpolation{T}}
        ie::H
        η2x::IT
    end
     
struct ConformalIEγν{T<:Real,  H <: IEγν{T}, IT <: AbstractInterpolation{T}}
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

function x_grid_ie(ie) # will this just work on IEγν? don't see why not...
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
    x_grid = x_grid_ie(ie) #Is this ever actually used? i.e. does x_grid_ie[1]=bg.x_grid[1]? looks like yes...
    xᵢ = first(x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=x_grid,
                dense=false,
                )
    return sol
end
function boltsolve(ie::IEν{T}, ode_alg=KenCarp4(); reltol=1e-6) where T 
    xᵢ = first(ie.bg.x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                # saveat=ie.bg.x_grid, 
                dense=false, #FIXME
                )
    return sol
end
function boltsolve(ie::IEallν{T}, ode_alg=KenCarp4(); reltol=1e-6) where T 
    xᵢ = first(ie.bg.x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false, #FIXME
                )
    return sol
end
function boltsolve(ie::IEγν{T}, ode_alg=KenCarp4(); reltol=1e-6) where T 
    xᵢ = first(ie.bg.x_grid)
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie!, u₀, (xᵢ , zero(T)), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false, #FIXME
                )
    return sol
end

function boltsolve_conformal(confie::ConformalIE{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    x_grid = η_grid_ie(ie,η2x,2048) #this is overkill/unoptomized but just to have something that decently agrees...
    xᵢ = first(x_grid) #to be consistent # again this does nothing...
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
    xᵢ = η2x( ie.bg.η[1] ) 
    Mpcfac = ie.bg.H₀*299792.458/100.
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie_conformal!, u₀, 
    (max(ie.bg.η[1]*Mpcfac,ie.bg.η(ie.bg.x_grid[1])*Mpcfac), 
    min(ie.bg.η[end]*Mpcfac,ie.bg.η(ie.bg.x_grid[end])*Mpcfac)),
        confie)
    sol = solve(prob, ode_alg, reltol=reltol,
                # saveat=ie.bg.η(x_grid)*Mpcfac,
                dense=false
                )
    return sol
end
function boltsolve_conformal(confie::ConformalIEallν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    xᵢ = ie.bg.x_grid[1]
    Mpcfac = ie.bg.H₀*299792.458/100.
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie_conformal!, u₀, 
        (max(ie.bg.η[1]*Mpcfac,ie.bg.η(ie.bg.x_grid[1])*Mpcfac), 
        min(ie.bg.η[end]*Mpcfac,ie.bg.η(ie.bg.x_grid[end])*Mpcfac)),
        confie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false
                )
    return sol
end
function boltsolve_conformal(confie::ConformalIEγν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    xᵢ = ie.bg.x_grid[1]#η2x( ie.bg.η[1] ) 
    # Mpcfac = ie.bg.H₀*299792.458/100.
    u₀ = initial_conditions(xᵢ, ie)
    prob = ODEProblem{true}(ie_conformal!, u₀, 
        # (max(ie.bg.η[1]*Mpcfac,ie.bg.η(ie.bg.x_grid[1])*Mpcfac), 
        # min(ie.bg.η[end]*Mpcfac,ie.bg.η(ie.bg.x_grid[end])*Mpcfac)),
        (max(ie.bg.η[1],ie.bg.η(ie.bg.x_grid[1])), 
        min(ie.bg.η[end],ie.bg.η(ie.bg.x_grid[end]))),
        confie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false
                )
    return sol
end

#FIXME: This copied code is unnecessary, but in the end we won't need 3 of these...
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
function ie_conformal!(du, u, confie::ConformalIEallν{T}, η) where T
    ie = confie.ie
    Mpcfac = ie.bg.H₀*299792.458/100.
    x = confie.η2x(η  / Mpcfac )
    ℋ = ie.bg.ℋ(x)
    ie!(du, u, ie, x)
    du .*= ℋ / Mpcfac  # account for dx/dη
    return nothing
end
function ie_conformal!(du, u, confie::ConformalIEγν{T}, η) where T
    ie = confie.ie
    Mpcfac = ie.bg.H₀*299792.458/100.
    # println("about to call η2x at ", η  / Mpcfac) # too many calls for this to be useful
    x = confie.η2x(η)#  / Mpcfac )
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
    ℓ_ν=0
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    # 𝒩 = OffsetVector(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1)) , 0:ℓ_ν)  # indexed 0 through ℓ_ν
    𝒩 = view(u, (2(ℓᵧ+1) + 1)) # only need dipole
    ℳ = OffsetVector(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq )) , 0:(ℓ_mν+1)*nq -1)  # indexed 0 through ℓ_mν
    Φ, δ, v, δ_b, v_b = view(u, ((2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5)) #getting a little messy...
    return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
end
function unpack(u, ie::IEallν{T, BasicNewtonian}) where T
    ℓᵧ =  ie.ℓ_γ
    nq = ie.nq
    ℓ_ν=0 #2 
    ℓ_mν = 0 #2 
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    # 𝒩 = OffsetVector(view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1)) , 0:ℓ_ν)  # indexed 0 through ℓ_ν
    𝒩 = view(u, (2(ℓᵧ+1) + 1)) # only need dipole
    # ℳ = OffsetVector(view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq )) , 0:(ℓ_mν+1)*nq -1)  # indexed 0 through ℓ_mν
    ℳ = view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq ))   # only need dipole (at all q)
    Φ, δ, v, δ_b, v_b = view(u, ((2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1 :(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5)) #getting a little messy...
    return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
end
function unpack(u, ie::IEγν{T, BasicNewtonian}) where T
    ℓᵧ = 1 #only monopole and dipole for both scalar temp and polzn
    nq = ie.nq
    ℓ_ν=0 #2 
    ℓ_mν = 0 #2 
    #NB we will keep the 0-1 indexing for Θ since it corresponds to multipoles, but massive neutrinos get 1->nq
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indexed 0 through ℓᵧ
    𝒩 = view(u, (2(ℓᵧ+1) + 1)) # only need dipole
    ℳ = view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq ))   # only need dipole (at all q)
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
function ie_unpack(u, ie::IEγν{T, BasicNewtonian}) where T
    ℓ_ν =  0
    ℓ_mν = 0
    nq = ie.nq
    Nᵧ = ie.Nᵧ₁+ie.Nᵧ₂+ie.Nᵧ₃ 
    ℓᵧ=1 #only monopole and dipole for both scalar temp and polzn
    Θ = OffsetArray(view(u, 1:(ℓᵧ+1),:), 0:ℓᵧ, 1:Nᵧ)  # indexed 0 through ℓᵧ, 1 through Nᵧ
    Θᵖ = OffsetArray(view(u, (ℓᵧ+2):(2ℓᵧ+2),:), 0:ℓᵧ, 1:Nᵧ)  # indexed 0 through ℓᵧ
    𝒩 = view(u, (2(ℓᵧ+1) + 1):(2(ℓᵧ+1)+ℓ_ν+1),:)   # indexed 0 through ℓ_ν
    ℳ = view(u, (2(ℓᵧ+1)+(ℓ_ν+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq ),:)   # indexed 0 through ℓ_mν
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
    Mpcfac = ie.bg.H₀*299792.458/100. 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x)/Mpcfac, ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)

    ℓ_ν = ie.ℓ_ν
    ℓ_mν =  ie.ℓ_mν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, ie)  # will be sweetened by .. syntax in 1.6
    # Θ[2] = ie.sΘ2(x)# call the spline, update Θ₂ at top since we do not evolve it
    Θ₂ = ie.sΘ2(x)# call the spline, update Θ₂ at top since we do not evolve it

    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ₂+#Θ[2]+
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
        # Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
        Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ₂ + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
        
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
    k, ℓ_ν, par, bg, ih, nq = ie.k, 0, ie.par, ie.bg, ie.ih, ie.nq
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
    Mpcfac = ie.bg.H₀*299792.458/100.
    𝒩₀ = ie.s𝒩₀(x)
    𝒩₂ = ie.s𝒩₂(x)

    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par)

    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩₂ #𝒩[2]#add rel quadrupole
                                  + σℳ / bg.ρ_crit /4
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩₀ #𝒩[0] #add rel monopole on this line
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
    𝒩′ = k/(3ℋₓ) * 𝒩₀ - 2*k/(3ℋₓ) *𝒩₂ + k/(3ℋₓ) *Ψ

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

function ie!(du, u, ie::IEallν{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓ_ν, ℓ_mν, par, bg, ih, nq = ie.k, 0, 0, ie.par, ie.bg, ie.ih, ie.nq #zeros here used to be 2s
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
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, ie)  # will be sweetened by .. syntax in 1.6
    
    
    #update pert vectors with splines
    𝒩₀ = ie.s𝒳₀[1](x) 
    𝒩₂ = ie.s𝒳₂[1](x)
    # WARNING no longer an offset array!
    ℳ₀ = zeros(T,nq)
    ℳ₂ = zeros(T,nq)
    for idx_q in 1:nq
        ℳ₀[idx_q] = ie.s𝒳₀[idx_q+1](x)
        ℳ₂[idx_q] = ie.s𝒳₂[idx_q+1](x)
    end
    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  @views ρ_σ(ℳ₀, ℳ₂, bg, a, par)

    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩₂    #𝒩[2]#add rel quadrupole
                                  + σℳ / bg.ρ_crit /4
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩₀ #𝒩[0] #add rel monopole on this line
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
        ℳ′[i_q+1] = k / (3ℋₓ) * ( q/ϵ * (ℳ₀[i_q+1] - 2ℳ₂[i_q+1])  - ϵ/q * Ψ  * df0)
    end

    # relativistic neutrinos (massless)
    𝒩′ = k/(3ℋₓ) * 𝒩₀ - 2*k/(3ℋₓ) *𝒩₂ + k/(3ℋₓ) *Ψ


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


function ie!(du, u, ie::IEγν{T, BasicNewtonian}, x) where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, ℓ_ν, ℓ_mν, par, bg, ih, nq = ie.k, 1, 0, 0, ie.par, ie.bg, ie.ih, ie.nq #zeros here used to be 2s
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    Mpcfac = ie.bg.H₀*299792.458/100. 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x)/Mpcfac, ih.τ′(x), ih.τ′′(x)
    #FIXME drop the unused things
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    csb² = ih.csb²(x)

    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  
    Θ′, Θᵖ′, 𝒩′, ℳ′, _, _, _, _, _ = unpack(du, ie)  
    
    #get perts from interpolators
    Θ₂ = ie.sΘ2(x)
    Π = ie.sΠ(x)
    𝒩₀ = ie.s𝒳₀[1](x) 
    𝒩₂ = ie.s𝒳₂[1](x)
    # WARNING no longer an offset array!
    ℳ₀ = zeros(T,nq)
    ℳ₂ = zeros(T,nq)
    for idx_q in 1:nq
        ℳ₀[idx_q] = ie.s𝒳₀[idx_q+1](x)
        ℳ₂[idx_q] = ie.s𝒳₂[idx_q+1](x)
    end
    #do the q integrals for massive neutrino perts (monopole and quadrupole)
    ρℳ, σℳ  =  @views ρ_σ(ℳ₀, ℳ₂, bg, a, par)

    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ₂
                                  + Ω_ν * 𝒩₂
                                  + σℳ / bg.ρ_crit /4
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩₀ 
        + a^(-2) * ρℳ / bg.ρ_crit
        )

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)

    # neutrinos (massive dipole, MB 57)
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
        ϵ,df0 = √(q^2 + (a*m_ν)^2), dlnf0dlnq(q,par)
        ℳ′[i_q+1] = k / (3ℋₓ) * ( q/ϵ * (ℳ₀[i_q+1] - 2ℳ₂[i_q+1])  - ϵ/q * Ψ  * df0)
    end

    # relativistic neutrinos (massless dipole)
    𝒩′ = k/(3ℋₓ) * 𝒩₀ - 2*k/(3ℋₓ) *𝒩₂ + k/(3ℋₓ) *Ψ

    # photons 
    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ₂ + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)

    # polarized photons
    Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    Θᵖ₂ = Π - Θᵖ[0] - Θ₂ #could drop this line but it makes things clearer
    Θᵖ′[1] = k / (3ℋₓ) * Θᵖ[0] - 2k / (3ℋₓ) * Θᵖ₂ + τₓ′ * Θᵖ[1] 

    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end



#FIXME: we don't actually need ANY of these if we are going to just start from a hierarchy call?
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
    k, ℓ_ν, par, bg, ih, nq = ie.k, 0, ie.par, ie.bg, ie.ih, ie.nq
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
    # Θ₀ = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    # Θ₂  = -8k / (15ℋₓ * τₓ′) * Θ[1]

    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    # Θᵖ[0] = (5/4) * Θ₂
    # Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ₂
    # Θᵖ[2] = (1/4) * Θ₂
    for ℓ in 3:ℓᵧ
        Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end
    # Θ[3] = -3/7 * k/(ℋₓ * τₓ′) * Θ₂
    # Θᵖ[3] = -3/7 * k/(ℋₓ * τₓ′) * Θᵖ[2]
    # for ℓ in 4:ℓᵧ
    #     Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
    #     Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    # end

    δ = 3/4 *(4Θ[0]) #the 4 converts δγ_MB -> Dodelson convention
    # δ = 3/4 *(4Θ₀) 
    δ_b = δ
    #we have that Θc = Θb = Θγ = Θν, but need to convert Θ = - k v (i absorbed in v)
    v = -3k*Θ[1]
    v_b = v

    # neutrino hierarchy
    # we need xᵢ to be before neutrinos decouple, as always
    # 𝒩[0] = Θ[0]
    𝒩₀ = Θ[0]
    # 𝒩₀ = Θ₀
    # 𝒩[1] = Θ[1]
    𝒩 = Θ[1]
    # 𝒩[2] = - (k^2 *ηₓ^2)/15 * 1 / (1 + 2/5 *f_ν) * Φ  / 2 #MB
    𝒩₂ = - (k^2 *ηₓ^2)/15 * 1 / (1 + 2/5 *f_ν) * Φ  / 2 #MB
    #FIXME^put the C here for consistency
    # for ℓ in 3:ℓ_ν
    #     𝒩[ℓ] = k/((2ℓ+1)ℋₓ) * 𝒩[ℓ-1] #standard truncation
    # end

    #massive neutrino hierarchy
    #It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
        ϵ = √(q^2 + (aᵢ*par.Σm_ν)^2)
        df0 = dlnf0dlnq(q,par)
        ℳ[0* nq+i_q] = -𝒩₀  *df0
        ℳ[1* nq+i_q] = -ϵ/q * 𝒩 *df0
        ℳ[2* nq+i_q] = -𝒩₂  *df0  #drop quadratic+ terms in (ma/q) as in MB
        for ℓ in 3:ℓ_mν #same scheme for higher-ell as for relativistic
            ℳ[ℓ* nq+i_q] = q / ϵ * k/((2ℓ+1)ℋₓ) * ℳ[(ℓ-1)*nq+i_q] #approximation of Callin06 (72), but add q/ϵ - leaving as 0 makes no big difference
        end
    end

    u[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end
#FIXME we don't actually need this, because we never initialize with the truncated hierarchy anwyways (these days)
function initial_conditions(xᵢ, ie::IEallν{T, BasicNewtonian}) where T
    k, ℓ_ν,ℓ_mν, par, bg, ih, nq = ie.k, 0, 0,ie.par, ie.bg, ie.ih, ie.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    ℓᵧ = ie.ℓ_γ
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)
    ℋₓ, _, ηₓ, τₓ′, _ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ are mutable views (see unpack)

    aᵢ² = exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
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

    #FIXME drop these
    # 𝒩[0] = Θ[0]
    𝒩₀ = Θ[0]
    # 𝒩[1] = Θ[1]
    𝒩 = Θ[1]
    # 𝒩[2] = - (k^2 *ηₓ^2)/15 * 1 / (1 + 2/5 *f_ν) * Φ  / 2 #MB
    𝒩₂ = - (k^2 *ηₓ^2)/15 * 1 / (1 + 2/5 *f_ν) * Φ  / 2 #MB

    #massive neutrino hierarchy
    #It is confusing to use Ψℓ bc Ψ is already the metric pert, so will use ℳ
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
        ϵ = √(q^2 + (aᵢ*par.Σm_ν)^2)
        df0 = dlnf0dlnq(q,par)
        # ℳ[0* nq+i_q] = -𝒩[0]  *df0
        # ℳ[1* nq+i_q] = -ϵ/q * 𝒩[1] *df0
        ℳ[1+i_q] = -ϵ/q * 𝒩 *df0
        # ℳ[2* nq+i_q] = -𝒩[2]  *df0  #drop quadratic+ terms in (ma/q) as in MB
    end

    u[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1:(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)] .= Φ, δ, v, δ_b, v_b  # write u with our variables
    return u
end

#FIXME we don't actually need this, because we never initialize with the truncated hierarchy anwyways (these days)
function initial_conditions(xᵢ, ie::IEγν{T, BasicNewtonian}) where T
    k, ℓᵧ, ℓ_ν,ℓ_mν, par, bg, ih, nq = ie.k, 1, 0, 0,ie.par, ie.bg, ie.ih, ie.nq
    Tν =  (par.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    q_pts = xq2q.(bg.quad_pts,logqmin,logqmax)
    u = zeros(T, 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5)
    Mpcfac = ie.bg.H₀*299792.458/100. 
    ℋₓ, _, ηₓ, τₓ′, _ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ)/Mpcfac, ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ are mutable views (see unpack)

    aᵢ² = exp(xᵢ)^2
    aᵢ = sqrt(aᵢ²)
    f_ν = 1/(1 + 1/(7*(3/3)*par.N_ν/8 *(4/11)^(4/3)))

    # metric and matter perturbations
    Φ = 1.0
    #choosing Φ=1 forces the following value for C, the rest of the ICs follow
    C = -( (15 + 4f_ν)/(20 + 8f_ν) )

    #trailing (redundant) factors are for converting from MB to Dodelson convention for clarity
    Θ[0] = -40C/(15 + 4f_ν) / 4
    Θ[1] = 10C/(15 + 4f_ν) * (k^2 * ηₓ) / (3*k)
    Θ₂ = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ₂
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ₂


    δ = 3/4 *(4Θ[0]) #the 4 converts δγ_MB -> Dodelson convention
    δ_b = δ
    #we have that Θc = Θb = Θγ = Θν, but need to convert Θ = - k v (i absorbed in v)
    v = -3k*Θ[1]
    v_b = v

    # neutrino hierarchy
    𝒩 = Θ[1]
    for (i_q, q) in zip(Iterators.countfrom(0), q_pts)
        ϵ = √(q^2 + (aᵢ*par.Σm_ν)^2)
        df0 = dlnf0dlnq(q,par)
        ℳ[1+i_q] = -ϵ/q * 𝒩 *df0
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
    Mpcfac = bg.H₀*299792.458/100. #FIXME won't need this anymore when k in correct units
    τ′,η = ih.τ′,bg.η #all splines of x
    y = k*( η(x)-η(x′) )/Mpcfac #Bessel argument
    IΘ2 = ( Θ0 - Φ′/ (-τ′(x′))  )*j2(y) - ( v_b   - ( k/bg.ℋ(x′) )*Ψ / (-τ′(x′)) )*j2′(y)  - Π*R2(y) / 2 
    return IΘ2
end

function _IΠ(x, x′,k, Π, bg)
    η = bg.η #all splines of x
    Mpcfac = bg.H₀*299792.458/100. 
    y = k*( η(x)-η(x′) )/Mpcfac  #Bessel argument
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
# /// IC Free streaming ///
# Relevant Bessel functions (ℓ=0,1,2)
#ℓ=0
j0(x) = (x > 0.01) ? sin(x)/x : 1 - x^2 /6 + x^4 /120 - x^6 /5040
j0′(x) = -j1(x)
#ℓ=1
j1(x) =  (x > 0.01) ?  (sin(x) - x*cos(x))/x^2 : x/3 - x^3 /30 + x^5 /840
R1(x) =  (x > 0.01) ? j1(x) - 3j2(x)/x : 2x/15 - 2x^3 /105 + x^5 /1260
#ℓ=2
# j2(x) = (x > 0.01) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^3 : x^2 /15 - x^4 /210 + x^6 /7560
# j2′(x) = (x > 0.01) ? ( -x*(x^2 -9)*cos(x) + (4x^2 -9)*sin(x) ) / x^4 : 2x /15 - 2x^3 /105 + x^5 /1260
# j2′′(x) = (x > 0.2) ? ( x*(5x^2 -36)*cos(x) + (x^4 - 17x^2 +36)*sin(x) ) / x^5 : 2/15 - 2x^2 /35 + x^4 /252 - x^6 /8910
# R2(x) = (x > 0.2) ? -( j2(x) + 3j2′′(x) ) / 2 : -1/5 + 11x^2 /210 -x^4 /280 +17x^4 /166320
# The W coupling kernel (sum truncated at ℓ=2)
W00(x) = j0(x)
W01(x) = j1(x)
W02(x) = j2(x)
W21(x) = -R1(x)
W22(x) = -R2(x)
function Wsum(x,𝒳ᵢ₀,𝒳ᵢ₁,𝒳ᵢ₂)
    𝒳ₛ₀ = W00(x)*𝒳ᵢ₀ - 3W01(x)*𝒳ᵢ₁ + 5W02(x)*𝒳ᵢ₂  #ℓ=0 ( use the subscript ₛ for streaming, this is the "free-streaming" piece)
    𝒳ₛ₂ = W02(x)*𝒳ᵢ₀ - 3W21(x)*𝒳ᵢ₁ + 5W22(x)*𝒳ᵢ₂ #ℓ=2
    return 𝒳ₛ₀, 𝒳ₛ₂
end

function get_Φ′_Ψ(u,hierarchy::Hierarchy{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, par, bg, nq = hierarchy.k, hierarchy.par, hierarchy.bg,hierarchy.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
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
    return Φ′,Ψ
end

# Get the Φ' and Ψ (copy function in ie file) from hierarchy
function get_Φ′_Ψ(u,ie::IEν{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    𝒩[0] = ie.s𝒩₀(x)
    𝒩[2] = ie.s𝒩₂(x)#WHY DO WE NEED THIS HERE BUT NOT IN PHOTONS? AND NOT FOR MONO?
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
    return Φ′,Ψ
end
function get_Φ′_Ψ(u,ie::IEallν{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    𝒩₀ = ie.s𝒳₀[1](x)
    𝒩₂ = ie.s𝒳₂[1](x)
    ℳ₀ = zeros(T,nq)
    ℳ₂ = zeros(T,nq)
    for idx_q in 1:nq#0:(nq-1)
        # ℳ[0*nq+idx_q] = ie.s𝒳₀[idx_q+2](x)
        # ℳ[2*nq+idx_q] = ie.s𝒳₂[idx_q+2](x)
        ℳ₀[idx_q] = ie.s𝒳₀[idx_q+1](x)
        ℳ₂[idx_q] = ie.s𝒳₂[idx_q+1](x)
    end
    # ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    ρℳ, σℳ  =  @views ρ_σ(ℳ₀, ℳ₂, bg, a, par)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩₂
                                  + σℳ / bg.ρ_crit /4
                                  )
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩₀
        + a^(-2) * ρℳ / bg.ρ_crit
        )
    return Φ′,Ψ
end

function fft_funcs(x, y, Φ′,Ψ, k,ℋ,q,m,𝕡)
    ϵ = (q^2 .+ exp.(2x)*m^2 ).^(1/2) #sqrt syntax doesn't work w/ bcast but put it back when undo bcast...
    q̃ = ϵ/q #convenience notation
    G₀ = ℋ .* q̃/k .* Φ′ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #for integrating in y #
    G₁ = -q̃.^2 .* Ψ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #
    K₀₀ = j0.(y) #1st index is ℓ 2nd index is derivative order
    K₀₁ = j0′.(y)
    K₂₀ = j2.(y) #
    K₂₁ = j2′.(y) #
    return G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁
end

function fft_integral(x, y,Φ′,Ψ,k,ℋ,q,m,𝕡,M) # for massive or massless neutrinos (𝒳=𝒩,ℳ)
    dy = y[2]-y[1]
    #  all ffts are performed in this function
    G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁ = fft_funcs(x,y, Φ′,Ψ, k,ℋ,q,m,𝕡) #
    # zero-pad the signals so convolution is not circular
    G₀,G₁ = [G₀; zeros(M-1)],[G₁; zeros(M-1)]
    K₀₀,K₀₁,K₂₀,K₂₁ = [K₀₀; zeros(M-1)],[K₀₁; zeros(M-1)],[K₂₀; zeros(M-1)],[K₂₁; zeros(M-1)] #
    # FFT the Gs, Ks
    G̃₀,G̃₁ = fft(G₀),fft(G₁)
    K̃₀₀, K̃₀₁, K̃₂₀, K̃₂₁ = fft(K₀₀),fft(K₀₁),fft(K₂₀),fft(K₂₁)#
    # Convolution theorem (iFFT pointwise product)
    𝒳₀ₓ = ifft(G̃₀.*K̃₀₀ .+ G̃₁.*K̃₀₁)[1:M]*dy 
    𝒳₂ₓ = ifft(G̃₀.*K̃₂₀ .+ G̃₁.*K̃₂₁)[1:M]*dy 
    return 𝒳₀ₓ,𝒳₂ₓ
end

function χ′z(a,q,m,tq_pts,tq_wts,𝕡)
    return q / (a * Bolt.ℋ_a(a,𝕡,tq_pts,tq_wts) * √(q^2 + (a*m)^2 ) )
end

function fft_ie(ie::IEν,perturb,M,m,q,i_q,u₀,x_grid)
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    # Set up the "neutrino horizon" and FFT abscissas
    # χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    χνs = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q1,𝕡.Σm_ν) for x in bg.x_grid]) #FIXME this won't work
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
    end
    _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
    if m==0 
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
    else
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
    end 
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)
end

function fft_ie_c(ie::IEν,perturb,M,m,q,i_q,u₀,x_grid) #FIXME add type decorators
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    # Set up the "neutrino horizon" and FFT abscissas
    # χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    χνs = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q1,𝕡.Σm_ν) for x in bg.x_grid])
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
    end

    _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie) 

    if m==0 
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
    else
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
    end 
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)#,
end


function fft_ie(ie::IEallν,perturb,M,u₀,x_grid)
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    #empty splines
    all_splines₀ = copy(ie.s𝒳₀)
    all_splines₂ = copy(ie.s𝒳₂)
    #explicitly do massless case
    χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),1.0,0.0,bg.quad_pts,bg.quad_wts,𝕡) for x in x_grid]) #bg.η
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources #FIXME this should probably happen outside of this function
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
    end

    𝒩₀ = ie.s𝒳₀[1](x_grid[1])
    _,_,𝒩₁,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒩₂ = ie.s𝒳₂[1](x_grid[1])
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩₁,𝒩₂)) #massless
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), 1.0,0.0,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    all_splines₀[1] = linear_interpolation(invx,𝒳₀)
    all_splines₂[1] = linear_interpolation(invx,𝒳₂)

    #massive case
    for i_q in 0:nq-1
        q = Bolt.xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),q,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts,𝕡) for x in x_grid])
        yyx = k.* (χνs .- χνs[1])
        dy=(yyx[end]-yyx[1])/(M-1)
        yy = yyx[1]:dy:yyx[end]
        invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
        
        # Get metric sources #FIXME this should probably happen outside of this function
        Φ′,Ψ = zeros(M),zeros(M)
        for j in 1:M
            Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
        end
        ℳ₀ = ie.s𝒳₀[2+i_q](x_grid[1])
        _,_,_, ℳ₁,_,_,_,_,_ =  unpack(u₀,ie)   
        ℳ₂ = ie.s𝒳₂[2+i_q](x_grid[1])

        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀,ℳ₁[1+i_q],ℳ₂)) #massive
        
        # Compute the new perts via FFT
        𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,𝕡.Σm_ν,𝕡,M)#,
        # Put it all together
        𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
        𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 

        all_splines₀[i_q+2] = linear_interpolation(invx,𝒳₀)
        all_splines₂[i_q+2] = linear_interpolation(invx,𝒳₂)

    end


    return invx, all_splines₀, all_splines₂
end

function fft_ie_c(ie::IEallν,perturb,M,u₀,x_grid) #FIXME add type decorators
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    all_splines₀ = copy(ie.s𝒳₀)
    all_splines₂ = copy(ie.s𝒳₂)
    #explicitly do massless case
    χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),1.0,0.0,bg.quad_pts,bg.quad_wts,𝕡) for x in x_grid])
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources #FIXME this should probably happen outside of this function
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
    end
    
    𝒩₀ = ie.s𝒳₀[1](x_grid[1])
    _,_,𝒩₁,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒩₂ = ie.s𝒳₂[1](x_grid[1])
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩₁,𝒩₂)) #massless
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), 1.0,0.0,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    all_splines₀[1] = linear_interpolation(invx,𝒳₀)
    all_splines₂[1] = linear_interpolation(invx,𝒳₂)

    #massive case
    for i_q in 0:nq-1
        q = Bolt.xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),q,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts,𝕡) for x in x_grid])
        yyx = k.* (χνs .- χνs[1])
        dy=(yyx[end]-yyx[1])/(M-1)
        yy = yyx[1]:dy:yyx[end]
        invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
        
        # Get metric sources #FIXME this should probably happen outside of this function
        Φ′,Ψ = zeros(M),zeros(M)
        for j in 1:M
            Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
        end
        ℳ₀ = ie.s𝒳₀[2+i_q](x_grid[1])
        _,_,_, ℳ₁,_,_,_,_,_ =  unpack(u₀,ie)   
        ℳ₂ = ie.s𝒳₂[2+i_q](x_grid[1])
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀,ℳ₁[1+i_q],ℳ₂)) #massive
        
        # Compute the new perts via FFT
        𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,𝕡.Σm_ν,𝕡,M)#,
        # Put it all together
        𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
        𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 

        all_splines₀[i_q+2] = linear_interpolation(invx,𝒳₀)
        all_splines₂[i_q+2] = linear_interpolation(invx,𝒳₂)

    end


    return invx, all_splines₀, all_splines₂
end

function h_boltsolve_flex(hierarchy::Hierarchy{T},  x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T
    prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (x_ini , x_fin), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false,
                )
    return sol
end

function boltsolve_flex(ie::IEν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
    prob = ODEProblem{true}(Bolt.ie!, u₀, (x_ini , x_fin), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false, 
                )
    return sol
end

function boltsolve_flex(ie::IEallν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
    prob = ODEProblem{true}(Bolt.ie!, u₀, (x_ini , x_fin), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false, 
                )
    return sol
end
function boltsolve_flex(ie::IEγν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
    prob = ODEProblem{true}(Bolt.ie!, u₀, (x_ini , x_fin), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false, 
                )
    return sol
end

function h_boltsolve_conformal_flex(confhierarchy::ConformalHierarchy{T},#FIXME we do't need this? {Hierarchy{T},AbstractInterpolation{T}},
    η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
    hierarchy = confhierarchy.hierarchy
    # Mpcfac = hierarchy.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(Bolt.hierarchy_conformal!, u₀, 
                            # (η_ini*Mpcfac , η_fin*Mpcfac),
                            (η_ini , η_fin),
                            confhierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,
    dense=false
    )
    return sol
end

function boltsolve_conformal_flex(confie::ConformalIEν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    # Mpcfac = ie.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(Bolt.ie_conformal!, u₀, 
                            # (η_ini*Mpcfac, η_fin*Mpcfac),
                            (η_ini, η_fin),
                            confie)
    sol = solve(prob, ode_alg, reltol=reltol,
    dense=false
    )
    return sol
end
function boltsolve_conformal_flex(confie::ConformalIEallν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    Mpcfac = ie.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(Bolt.ie_conformal!, u₀, 
                            (η_ini*Mpcfac, η_fin*Mpcfac),
                            confie)
    sol = solve(prob, ode_alg, reltol=reltol,
    dense=false
    )
    return sol
end
function boltsolve_conformal_flex(confie::ConformalIEγν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    # Mpcfac = ie.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(Bolt.ie_conformal!, u₀, 
                            # (η_ini*Mpcfac, η_fin*Mpcfac),
                            (η_ini, η_fin),
                            confie)
    sol = solve(prob, ode_alg, reltol=reltol,
    dense=false
    )
    return sol
end

function iterate_fft(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
    M, reltol,x_ini, x_fin,u0,m,q,i_q) where T
    # 𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
    ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒩₀_km1, 𝒩₂_km1,
                    ℓᵧ, ℓ_mν, n_q)
    perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie(ie_k_late,perturb_k_late,M,m,q,i_q,
                        u0,perturb_k_late.t) 
    return xx,𝒩₀_k,𝒩₂_k,perturb_k_late
end

function iterate_fft_c(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
    M, reltol,η_ini, η_fin,u0,m,q,i_q) where T
    # 𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
    ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒩₀_km1, 𝒩₂_km1,
                    ℓᵧ, ℓ_mν, n_q)
    ie_k_conf_late_c = ConformalIEν(ie_k_late,η2x);
    perturb_k_late_c = boltsolve_conformal_flex(ie_k_conf_late_c, η_ini, η_fin, u0; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie_c(ie_k_conf_late_c.ie,perturb_k_late_c,M,m,q,i_q,
                        u0,
                        # η2x(perturb_k_late_c.t/Mpcfac)
                        η2x(perturb_k_late_c.t)
                        ) 
    return xx,𝒩₀_k,𝒩₂_k,perturb_k_late_c
end

function iterate_fft_allν(𝒳₀_km1,𝒳₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, n_q,
    M, reltol,x_ini, x_fin,u0) where T
    ie_k_late = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
                     𝒳₀_km1,𝒳₂_km1,
                    ℓᵧ, n_q)
    perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
    xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k_late,perturb_k_late,M,
                        u0,perturb_k_late.t) 
    return xx,𝒳₀_k,𝒳₂_k,perturb_k_late
end

function iterate_fft_allν_c(𝒳₀_km1,𝒳₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, n_q,
    M, reltol,η_ini, η_fin,u0) where T
    ie_k_late = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒳₀_km1,𝒳₂_km1,
                    ℓᵧ, n_q)
    ie_k_conf_late_c = ConformalIEallν(ie_k_late,η2x);
    perturb_k_late_c = boltsolve_conformal_flex(ie_k_conf_late_c, η_ini, η_fin, u0; reltol=reltol)
    xx,𝒳₀_k,𝒳₂_k = fft_ie_c(ie_k_conf_late_c.ie,perturb_k_late_c,M,
                        u0,
                        # η2x(perturb_k_late_c.t/Mpcfac)
                        η2x(perturb_k_late_c.t)
                        ) 
    return xx,𝒳₀_k,𝒳₂_k,perturb_k_late_c
end

#---------------------------------#
# Itersolves
#---------------------------------#
function itersolve_fft(Nₖ::Int,ie_0::IEν{T},M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
    𝒩₀_0,𝒩₂_0 = ie_0.s𝒩₀,ie_0.s𝒩₂
    𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
    perturb_k = nothing
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft(𝒩₀_k,𝒩₂_k,ie_0.par,ie_0.bg,ie_0.ih,
                                   ie_0.k,ie_0.ℓ_γ,ie_0.ℓ_mν,ie_0.nq,
                                   M,reltol,x_ini,x_fin,u0,m,q,i_q)
    end
    return xx_k, 𝒩₀_k,𝒩₂_k,perturb_k
end
#ctime version
function itersolve_fft(Nₖ::Int,ie_0_c::ConformalIEν{T},M::Int,η_ini, η_fin,u0;reltol=1e-6) where T
    𝒩₀_0,𝒩₂_0 = ie_0_c.ie.s𝒩₀,ie_0_c.ie.s𝒩₂
    𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
    perturb_k = nothing
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft_c(𝒩₀_k,𝒩₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.ℓ_mν,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0,m,q,i_q)
    end
    return xx_k,𝒩₀_k,𝒩₂_k,perturb_k
end

function itersolve_fft(Nₖ::Int,ie_0::IEallν{T},M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
    𝒳₀_0,𝒳₂_0 = ie_0.s𝒳₀,ie_0.s𝒳₂
    𝒳₀_k,𝒳₂_k = 𝒳₀_0,𝒳₂_0 
    perturb_k = nothing
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate_fft_allν(𝒳₀_k,𝒳₂_k,ie_0.par,ie_0.bg,ie_0.ih,
                                   ie_0.k,ie_0.ℓ_γ,ie_0.nq,
                                   M,reltol,x_ini,x_fin,u0)
    end
    return xx_k,𝒳₀_k,𝒳₂_k,perturb_k
end
#ctime version
function itersolve_fft(Nₖ::Int,ie_0_c::ConformalIEallν{T},M::Int,η_ini, η_fin,u0;reltol=1e-6) where T
    𝒳₀_0,𝒳₂_0 = ie_0_c.ie.s𝒳₀,ie_0_c.ie.s𝒳₂
    𝒳₀_k,𝒳₂_k = 𝒳₀_0,𝒳₂_0
    perturb_k = nothing
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate_fft_allν_c(𝒳₀_k,𝒳₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0)
    end
    return xx_k,𝒳₀_k,𝒳₂_k,perturb_k
end


# Helper functon for switch
function get_switch_u0(η,hierarchy_conf) #Input is η of the switch
    # This function assumes truncated hierarchies for all neutrinos (but not yet photons)
    hierarchy,bg = hierarchy_conf.hierarchy,hierarchy_conf.hierarchy.bg
    Mpcfac = bg.H₀*299792.458/100.
    switch_idx = argmin(abs.(bg.η .-η)) #for now we use the bg to find the switch
    #solve the split ode
    ℓᵧ,ℓ_ν,n_q = hierarchy.ℓᵧ,hierarchy.ℓ_ν, hierarchy.nq
    pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
    # \/ we want to report this timing to get a full picture of total time (early+late)
    sol_early_c = h_boltsolve_conformal_flex(hierarchy_conf, bg.η[1], bg.η[switch_idx],  initial_conditions(bg.x_grid[1], hierarchy));
    
    # Get the new initial conditions
    u0_ie_c = zeros(2(ℓᵧ+1) + (0+1) + (0+1)*n_q + 5);
    # The first split will be the same
    for i in  1:2(ℓᵧ+1) #up to massless ν quadrupole (w/photon hierarchy #FIXME)
        u0_ie_c[i] = sol_early_c.u[end][i]
    end
    #set the massless neutrino dipole
    u0_ie_c[2(ℓᵧ+1)+1] = sol_early_c.u[end][2(ℓᵧ+1)+2]

    #This does nothing??
    # for i in  2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+1 #skip the higher massless hierarchy multipoles
    #     down_shift = i-(ℓ_ν-2) #shift the hierarchy index down by the number of multipoles we skip in the ie
    #     u0_ie_c[down_shift] = sol_early_c.u[end][i]
    # end

    #massive neutrinos, now we just do the dipole again

    #TODO Since this is contiguous, can combine it with the above loop
    # Do the same for massive neutrinos, which are arranged as [q1,q2,...,qnq]_ℓ=0, [q1,q2,...,qnq]_ℓ=1, ..., []_ℓ=ℓ_mν
    # for i_ℓ in 1:3 #we fill all the multipoles as we usually would up to the quadrupole
    #     for i in 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*(i_ℓ-1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+i_ℓ*n_q 
    #         down_shift = i-(ℓ_ν-2)
    #         u0_ie_c[down_shift] = sol_early_c.u[end][i]
    #     end
    # end

    # start at the dipole first q idx, go up to the last dipole q idx (in the hierarchy)   
    for i in 1:n_q 
        # down_shift = i-(ℓ_ν+1 - 1)-(ℓ_mν+1 - 1)
        u0_ie_c[2(ℓᵧ+1)+1+i] = sol_early_c.u[end][2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1+i]
    end

    for i in  1:5 #skip the higher massless hierarchy multipoles
        # down_shift = i-(ℓ_ν+1-1)-n_q*(ℓ_mν+1-1) #shift down by the number of multipoles we skip
        u0_ie_c[2(ℓᵧ+1)+1+n_q+i] = sol_early_c.u[end][pertlen-5+i]
    end

    return u0_ie_c
end

#-------------------------------------------------------------------------------------------
#---------------------                  Final Iteration functions           ----------------
#-------------------------------------------------------------------------------------------

function x_grid_ie(ie,x_ini,x_fin) # We don't want this to use the background...
    bg,ih,k = ie.bg,ie.ih,ie.k
    # x_ini,x_fin = ...
    # assert(bg.η(ie.x_grid[1]) <= 2π/k) #check that we start earlier than horizon crossing
    #If we start before horizon entry we can still use the below scheme for setting x trapz points
    # and we should make Nᵧ₁ small
    #Actually, is this fine?
    # digits=3
    # Three phases: 
    # 1. Pre-horizon entry:
    Mpcfac = ie.bg.H₀*299792.458/100.
    xhor = bg.x_grid[argmin(abs.(k .* bg.η /Mpcfac .- 2π))] #horizon crossing ish
    # x_ph_i, x_ph_f, n_ph = bg.x_grid[1], xhor, ie.Nᵧ₁ #10.
    x_ph_i, x_ph_f, n_ph = x_ini, xhor, ie.Nᵧ₁ 
    dx_ph = (x_ph_f-x_ph_i)/(n_ph-1)
    # x_ph = -20.:dx_ph:x_ph_f
    # x_ph = -20. .+ dx_ph*collect(0:1:n_ph-1)
    # x_ph = x_ph_i .+ dx_ph*collect(0:1:n_ph-1)
    x_ph = collect(range(x_ph_i,x_ph_f,length=n_ph-1))

    # 2. Wiggly time (recomb):
    xdec = bg.x_grid[argmin(abs.( -ih.τ′ .* bg.ℋ .*bg.η /Mpcfac .- 1))] #decoupling ish
    x_rc_f, n_rc = xdec, ie.Nᵧ₂ #100
    dx_rc = (x_rc_f-x_ph_f)/n_rc
    # x_rc = x_ph_f .+ dx_rc * collect(1:1:n_rc)
    x_rc = range(x_ph_f,xdec,length=n_rc-1)

    # 3. Post-recomb:
    n_pr = ie.Nᵧ₃ #50
    # dx_pr = (bg.x_grid[end] -x_rc_f)/n_pr
    dx_pr = (x_fin -x_rc_f)/n_pr
    # x_pr = x_rc_f .+ dx_pr* collect(1:1:n_pr )
    x_pr = range(x_ph_f,xdec,length=n_rc-1)
    # println("1,2,3 = ($(x_ph_i),$(x_ph_f),$(x_rc_f))")
    # println("dx1,dx2,dx3 = ($(dx_ph),$(dx_rc),$(dx_pr))")
    x_sparse = vcat(x_ph,x_rc,x_pr)


    x_sparse = unique([range(x_ph_i,x_ph_f,n_ph+1);
                       range(x_ph_f,x_rc_f,n_rc+1); 
                       range(x_rc_f,x_fin,n_pr)])

    return x_sparse
end


function η_grid_ie(ie,η_ini,η_fin,η2x) # We don't want this to use the background...
    bg,ih,k = ie.bg,ie.ih,ie.k
    η_grid = range(η_ini,η_fin,1024)
    τ′ = ih.τ′(η2x(η_grid))
    ℋ = bg.ℋ(η2x(η_grid))
    # Three phases: 
    # 1. Pre-horizon entry:
    Mpcfac = ie.bg.H₀*299792.458/100.
    ηhor = η_grid[argmin(abs.(k .* η_grid /Mpcfac .- 2π))] #horizon crossing ish
    # x_ph_i, x_ph_f, n_ph = bg.x_grid[1], xhor, ie.Nᵧ₁ #10.
    η_ph_i, η_ph_f, n_ph = η_ini, ηhor, ie.Nᵧ₁ 
    # 2. Wiggly time (recomb):
    ηdec = η_grid[argmin(abs.( -τ′ .* ℋ .*η_grid /Mpcfac .- 1))] #decoupling ish
    η_rc_f, n_rc = ηdec, ie.Nᵧ₂ 
    # 3. Post-recomb:
    n_pr = ie.Nᵧ₃ 
    η_sparse = unique([range(η_ph_i,η_ph_f,n_ph+1);
                       range(η_ph_f,η_rc_f,n_rc+1); 
                       range(η_rc_f,η_fin,n_pr)])

    return η_sparse
end


# This should just work...
function g_weight_trapz_ie(i,x_grid,ie::IEγν{T},Φ′,Ψ,Θ₀,Π,v_b) where T
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
    Θ2ᵢ += (_IΘ2(xᵢ,x_grid[i],k,0.0,Θ₀[i],v_b[i],Φ′[i],Ψ[i],ie.ih,ie.bg)*Wⱼ⁺
           + _IΘ2(xᵢ,x_grid[i-1],k,Π[i-1],Θ₀[i-1],v_b[i-1],Φ′[i-1],Ψ[i-1],ie.ih,ie.bg)*Wⱼ)
    Πᵢ += _IΠ(xᵢ,x_grid[i-1],k,Π[i-1],ie.bg)*Wⱼ
    #Kamionkowski integration scheme for handling x′ = x at each x (the implicit timestep)
    Πᵢ = (Πᵢ + Θ2ᵢ) / ( 1 - 7/10 * Wⱼ⁺)
    Θ2ᵢ = Θ2ᵢ + ( i <length(x_grid) ? Π[i+1] : 0. )/10 * Wⱼ⁺ #if i+1>length(x_grid), return 0 for oob array
    return Θ2ᵢ,Πᵢ
end

function fft_ie(ie::IEγν,M,u₀,x_grid,sΦ′,sΨ)
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    #empty interpolators
    all_splines₀ = copy(ie.s𝒳₀)
    all_splines₂ = copy(ie.s𝒳₂)

    # Do the massless case
    χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),1.0,0.0,bg.quad_pts,bg.quad_wts,𝕡) for x in x_grid]) #bg.η
    yyx = k.* (χνs .- χνs[1])
    # println("x_grid[1]", x_grid[1])
    # println("x_grid[end]", x_grid[end])
    # println("χνs[end]", χνs[end])
    # println("k", k)
    # println("yyx[end]: ", yyx)
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources #FIXME this should probably happen outside of this function
    #^It can only happen outside if we save interpolators for the metric and call them,
    #which is probably what we should be doing...
    # Φ′,Ψ = zeros(M),zeros(M)
    # for j in 1:M
    #     Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
    # end
    # println("invx[1]", invx[1])
    Φ′,Ψ = sΦ′(invx),sΨ(invx)

    # Free-streaming piece
    𝒩₀ = ie.s𝒳₀[1](x_grid[1])
    _,_,𝒩₁,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒩₂ = ie.s𝒳₂[1](x_grid[1])
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩₁,𝒩₂)) #massless
    # FFT piece
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), 1.0,0.0,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    all_splines₀[1] = linear_interpolation(invx,𝒳₀)
    all_splines₂[1] = linear_interpolation(invx,𝒳₂)


    # Loop over q for the massive case
    for i_q in 0:nq-1
        # Neutrino horizon 
        #FIXME neutrino horizon function
        q = Bolt.xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),q,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts,𝕡) for x in x_grid])
        yyx = k.* (χνs .- χνs[1])
        dy=(yyx[end]-yyx[1])/(M-1)
        yy = yyx[1]:dy:yyx[end]
        invx = linear_interpolation(yyx,x_grid).(yy)  #NB, mass and q depenent
        
        # Get metric sources #FIXME this should probably happen outside of this function
        # Φ′,Ψ = zeros(M),zeros(M)
        # for j in 1:M
        #     Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
        # end
        Φ′,Ψ = sΦ′(invx),sΨ(invx)

        # Free-streaming piece
        ℳ₀ = ie.s𝒳₀[2+i_q](x_grid[1])
        _,_,_, ℳ₁,_,_,_,_,_ =  unpack(u₀,ie)   
        ℳ₂ = ie.s𝒳₂[2+i_q](x_grid[1])
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀,ℳ₁[1+i_q],ℳ₂)) 
        
        # FFT piece
        𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,𝕡.Σm_ν,𝕡,M)#,

        # Put it all together
        𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
        𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 

        # Update massive neutrinos
        all_splines₀[i_q+2] = linear_interpolation(invx,𝒳₀)
        all_splines₂[i_q+2] = linear_interpolation(invx,𝒳₂)

    end


    return invx, all_splines₀, all_splines₂
end

