using CUDA,Adapt,Interpolations
using DiffEqGPU, OrdinaryDiffEq
using Pkg
Pkg.activate("/pscratch/sd/j/jsull/julia/Bolt.jl")
using Bolt
using ForwardDiff
#using Plots
using BenchmarkTools
using Setfield,StaticArrays,OffsetArrays

# bg/ion setup
𝕡 = CosmoParams()
const n_q=15
const nq=15

logqmin,logqmax = -6,-1
bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
kmin,kmax= 0.1bg.H₀*100,5000bg.H₀
k_grid = log10_k(kmin,kmax,33)

const ℓᵧ=25
const ℓ_ν=25
const ℓ_mν=20
reltol=1e-5 #cheaper  rtol
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(bg.x_grid))
ℳρ,ℳσ = zeros(length(bg.x_grid)),zeros(length(bg.x_grid)) #arrays for the massive neutrino integrated perts
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k_grid[1], ℓᵧ, ℓ_ν, ℓ_mν,n_q)
# storage = CUDA.CuArrayAdaptor{Mem.DeviceBuffer}() # for Float32 CuArray
# hierarchy = adapt(storage, Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k_grid[1], 25,25,20,15))

#function k_boltsolve(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6) where T
xᵢ = first(hierarchy.bg.x_grid)
u₀ = cu(Bolt.initial_conditions(xᵢ, hierarchy))
du = zero(u₀)


function Adapt.adapt_structure(to, itp::Interpolations.BSplineInterpolation{T,N,<:Any,IT,Axs}) where {T,N,IT,Axs}
    coefs = Adapt.adapt_structure(to, itp.coefs)
    Tcoefs = typeof(coefs)
    Interpolations.BSplineInterpolation{T,N,Tcoefs,IT,Axs}(coefs, itp.parentaxes, itp.it)
end

function Adapt.adapt_structure(to, itp::Interpolations.ScaledInterpolation{T,N,ITPT,IT,<:Any}) where {T,N,ITPT,IT}
    s = Adapt.adapt_structure(to,itp.itp)
    Titp = typeof(s)
    ranges = Adapt.adapt_structure(to,itp.ranges)
    RT=typeof(ranges)
    Interpolations.ScaledInterpolation{T,N,Titp,IT,RT}(s,ranges)
end

function gpu_unpack(u)  #use Marius' trick for the ntuples to avoid size limits on tuples
    Θ = OffsetVector(SVector(ntuple(i -> u[i], Val(ℓᵧ+1))), 0:ℓᵧ)
    Θᵖ = OffsetVector(SVector(ntuple(i -> u[i+(ℓᵧ+1)], Val(ℓᵧ+1))), 0:ℓᵧ)
    𝒩 = OffsetVector(SVector(ntuple(i -> u[i+2(ℓᵧ+1)], Val(ℓ_ν+1))), 0:ℓ_ν)
    ℳ  = OffsetVector(SVector(ntuple(i -> u[i+2(ℓᵧ+1)+(ℓ_ν+1)], Val((ℓ_mν+1)*nq))), 0:(ℓ_mν+1)*nq -1)
    Φ, δ, v, δ_b, v_b = SVector(ntuple(i -> u[i+2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq], Val(5)))
    return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
end

##

#adapts
Adapt.@adapt_structure Background
Adapt.@adapt_structure IonizationHistory
Adapt.@adapt_structure Hierarchy #give this a shot not sure it's gonna work


function f_kernel!(du, u, h, x)
    # get all the data
    nq = n_q
    k, par,bg,ih = h.k, h.par, h.bg,h.ih
    Ω_r, Ω_b, Ω_m, N_ν, m_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, par.Σm_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′, csb² = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x),ih.csb²(x)
    a = x2a(x)
    Tν =  (N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(par) *Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    
    R = 4Ω_r / (3Ω_b * a)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    # ρℳ, σℳ  =  Bolt.ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    ϵx(x, am) = √(Bolt.xq2q(x,logqmin,logqmax)^2 + (am)^2)
    Iρ(x) = Bolt.xq2q(x,logqmin,logqmax)^2  * ϵx(x, a*m_ν) * Bolt.f0(Bolt.xq2q(x,logqmin,logqmax),par) / Bolt.dxdq(Bolt.xq2q(x,logqmin,logqmax),logqmin,logqmax)
    Iσ(x) = Bolt.xq2q(x,logqmin,logqmax)^2  * (Bolt.xq2q(x,logqmin,logqmax)^2 /ϵx(x, a*m_ν)) * f0(Bolt.xq2q(x,logqmin,logqmax),par) / Bolt.dxdq(Bolt.xq2q(x,logqmin,logqmax),logqmin,logqmax)
    xq,wq = bg.quad_pts,bg.quad_wts


    # do the unpack
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = gpu_unpack(u)
    Θ′, Θᵖ′, 𝒩′, ℳ′, Φ′, δ′, v′, δ_b′, v_b′ = gpu_unpack(du)

    ρℳ, σℳ  =  0.,0.
    for i in 1:nq #have to un-broadcast this...
        ρℳ += 4π*Iρ(xq[1])*ℳ[0*nq+i-1]*wq[i]
        σℳ += 4π*Iσ(xq[i])*ℳ[2*nq+i-1]*wq[i]
    end

    #start setting the perturbations
    # metric
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
    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)

    # relativistic neutrinos (massless)
    @set! 𝒩′[0] = -k / ℋₓ * 𝒩[1] - Φ′ #for some reason need set here...
    𝒩′[1] = k/(3ℋₓ) * 𝒩[0] - 2*k/(3ℋₓ) *𝒩[2] + k/(3ℋₓ) *Ψ
    for ℓ in 2:(ℓ_ν-1)
        𝒩′[ℓ] =  k / ((2ℓ+1) * ℋₓ) * ( ℓ*𝒩[ℓ-1] - (ℓ+1)*𝒩[ℓ+1] )
    end
    𝒩′[ℓ_ν] =  k / ℋₓ  * 𝒩[ℓ_ν-1] - (ℓ_ν+1)/(ℋₓ *ηₓ) *𝒩[ℓ_ν]

    # photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    @set! Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    for ℓ in 2:(ℓᵧ-1)
        Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * Bolt.δ_kron(ℓ, 2) / 10)
    end

    # # polarized photons
    @set! Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓᵧ-1)
        Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * Bolt.δ_kron(ℓ, 2) / 10)
    end

    # # photon boundary conditions: diffusion damping 
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[ℓᵧ]

    # massive neutrinos
    #FIXME ℳ′ assignment does not work in this loop for some reason??
    for i_q in 1:nq
        q = Bolt.xq2q(bg.quad_pts[i_q] ,logqmin,logqmax)
        ϵ = √(q^2 + (a*m_ν)^2)
        df0 = dlnf0dlnq(q,par)
        du[2(ℓᵧ+1)+(ℓ_ν+1) + 0*nq + i_q] = - k / ℋₓ *  q/ϵ * ℳ[1* 10+i_q-1]  + Φ′ * df0 #ℳ′[0*nq+i_q-1]
        du[2(ℓᵧ+1)+(ℓ_ν+1) + 1* nq+i_q] = k / (3ℋₓ) * ( q/ϵ * (ℳ[0* nq+i_q] - 2ℳ[2* nq+i_q])  - ϵ/q * Ψ  * df0) #ℳ′[1* nq+i_q]
        for ℓ in 2:(ℓ_mν-1)
            du[2(ℓᵧ+1)+(ℓ_ν+1) + ℓ* nq+i_q]=  k / ℋₓ * q / ((2ℓ+1)*ϵ) * ( ℓ*ℳ[(ℓ-1)* nq+i_q-1] - (ℓ+1)*ℳ[(ℓ+1)* nq+i_q-1] ) #ℳ′[ℓ* nq+i_q]
        end
        du[2(ℓᵧ+1)+(ℓ_ν+1) + ℓ_mν* nq+i_q] =  q / ϵ * k / ℋₓ * ℳ[(ℓ_mν-1)* nq+i_q-1] - (ℓ_mν+1)/(ℋₓ *ηₓ) *ℳ[(ℓ_mν)* nq+i_q-1] #ℳ′[ℓ_mν* nq+i_q]  MB (58) similar to rel case but w/ q/ϵ
    end

    # for some reason OffsetArray values are not mutating...do it by hand...
    for i in 1:(ℓᵧ+1)
        du[i] = Θ′[i-1]
        du[(ℓᵧ+1)+i] = Θᵖ′[i-1]
    end
    for i in 1:(ℓ_ν+1)
        du[2(ℓᵧ+1)+i] = 𝒩′[i-1]
    end
    # See above 
    # for i in 1:(ℓ_mν+1)
    #     for j in 1:nq
    #         du[2(ℓᵧ+1)+(ℓ_ν+1)+i*nq + j-1] = ℳ′[0* nq+i-1]
    #     end
    # end

    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+1] = Φ′
    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+2] = δ′
    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+3] = v′
    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+4] = δ_b′
    du[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq+5] = v_b′ 

   return nothing
end

#check these 
f_kernel!(du,u₀,cu(hierarchy),xᵢ) 
@cuda f_kernel!(du,u₀,cu(hierarchy),xᵢ)

#this used to be setting up ensemble problem but that didn't work...
function gpu_hierarchy!(du, u, hierarchy,x) 
    f_kernel!(du,u,cu(hierarchy),x) 
end

gpu_ensemble_hierarchy!(du,u₀,k_grid[1,:],xᵢ)
du
gpu_prob = ODEProblem{true}(gpu_hierarchy!, u₀, (xᵢ,zero(typeof(xᵢ))), hierarchy)

xᵢf = Float32(xᵢ)
gpu_prob = ODEProblem{true}(f_kernel!, u₀, (xᵢf,zero(typeof(xᵢf))),cu(hierarchy))
sol = solve(gpu_prob, KenCarp4(), reltol=1e-4,saveat=Float32.(hierarchy.bg.x_grid), dense=false,)
#^When I run this line I get a warning about setting scalars


# Try something simpler...Euler should have to work? But on first step says unstable...must be something wrong when calling the hierarachy function
xᵢ = first(hierarchy.bg.x_grid)
u₀ = cu(Bolt.initial_conditions(xᵢ, hierarchy))
du = zero(u₀)
sol = solve(gpu_prob, Euler(),dt=Float32(1e-3),saveat=Float32.(hierarchy.bg.x_grid), dense=false,)

#This has to work at some level? Why does it produce nans after a handful of steps?
N=3
for i in 1:N
    f_kernel!(du,u,hierarchy,x)
    println(x," du: ", du[1], " u: ", u[1] )
    println(u)
    @CUDA.sync u .= u .+ dt*du
    x=x+dt
    for j in length(u) 
        u[j] = u[j] + dt*du[j] 
        res[j,i] = u[j] + dt*du[j] 
    end
    u=u+dt*du
end
