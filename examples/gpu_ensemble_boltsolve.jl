using CUDA
using DiffEqGPU, OrdinaryDiffEq
using Pkg
Pkg.activate("/pscratch/sd/j/jsull/julia/Bolt.jl")
using Bolt
using ForwardDiff
#using Plots
using BenchmarkTools

# bg/ion setup
𝕡 = CosmoParams()
n_q=15
logqmin,logqmax = -6,-1
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
kmin,kmax= 0.1bg.H₀*100,5000bg.H₀
k_grid = log10_k(kmin,kmax,33)

ℓᵧ=50
ℓ_ν=50
ℓ_mν=20
reltol=1e-5 #cheaper  rtol
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(bg.x_grid))
ℳρ,ℳσ = zeros(length(bg.x_grid)),zeros(length(bg.x_grid)) #arrays for the massive neutrino integrated perts
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k_grid[1], ℓᵧ, ℓ_ν, ℓ_mν,n_q)

#function k_boltsolve(hierarchy::Hierarchy{T}, ode_alg=KenCarp4(); reltol=1e-6) where T
xᵢ = first(hierarchy.bg.x_grid)
u₀ = Bolt.initial_conditions(xᵢ, hierarchy)
T=Float32 #will not need this later when T is inferred from the function
prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
#to create the en
prob_func = (prob,i,repeat) -> remake(prob,p=Hierarchy(prob.p.integrator,
                                                       prob.p.par,
                                                       prob.p.bg,
                                                       prob.p.ih,
                                                       k_grid[i],
                                                       prob.p.ℓᵧ,
                                                       prob.p.ℓ_ν,
                                                       prob.p.ℓ_mν,
                                                       prob.p.nq)
                                     )
eprob = EnsembleProblem(prob, prob_func = prob_func, safetycopy=false)
#@time sol = solve(prob,KenCarp4(),EnsembleGPUArray(),trajectories=10,saveat=1.0f0) this does not work for some reason...
#^I am matching the syntax since the below works without trying to use the GPU
@time esol = solve(eprob,KenCarp4(),EnsembleThreads(),trajectories=length(k_grid),saveat=bg.x_grid,reltol=reltol,dense=false)
@time for i in length(k_grid)
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k_grid[i], ℓᵧ, ℓ_ν, ℓ_mν,n_q)
    prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    solve(prob,KenCarp4(),saveat=bg.x_grid,reltol=reltol,dense=false)(bg.x_grid)
end

#println(esol[end]/sols)
