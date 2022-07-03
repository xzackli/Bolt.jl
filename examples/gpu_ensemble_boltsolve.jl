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


# For the GPUEnsemble, wrap the hierarchy inside a function that takes k as an argument
# Uses global knowledge of hierarchy arguments for cosntruction
function gpu_ensemble_hierarchy!(du, u, k, x) #where T 
    #build the hierarchy
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
    #call the derivative function 
    Bolt.hierarchy!(du,u,hierarchy,x) 
end
#test this function on the first step
println("pre type of u0, ", typeof(u₀))
println("pre type of kgrid, ", typeof(k_grid))
println("pre type of xi, ", typeof(xᵢ))
u₀ = T.(u₀)
k_grid = T.(k_grid)
k_grid = reshape(k_grid,(size(k_grid)...,1)) #For some reason GPU very specifically wants array of 1 size arrays rather than array of floats...
println("post type of u0, ", typeof(u₀))
println("post type of kgrid, ", typeof(k_grid))
dutest = zero(u₀)
testderiv = gpu_ensemble_hierarchy!(dutest,u₀,k_grid[1],xᵢ)
println("Succesfully evaluted gpu derivative")
println("test deriv: ", testderiv)
println("one element",typeof(k_grid[1]))

gpu_prob = ODEProblem{true}(gpu_ensemble_hierarchy!, u₀, (xᵢ , zero(T)), k_grid[1,:])
gpu_prob_func = (gpu_prob,i,repeat) -> remake(gpu_prob,p=k_grid[i,:])
eprob = EnsembleProblem(gpu_prob,prob_func=gpu_prob_func,safetycopy=false)
#@time esol = solve(eprob,KenCarp4(),EnsembleThreads(),trajectories=length(k_grid),saveat=bg.x_grid,reltol=reltol,dense=false)
#@time esol = solve(eprob,Tsit5(),EnsembleThreads(),trajectories=length(k_grid),saveat=bg.x_grid,reltol=reltol,dense=false) #this runs but aborts the solve due to instability - maybe not high enough ellmax...
#println("solved thread ensemble")
#@time gsol = solve(eprob,KenCarp4(),EnsembleGPUArray(),trajectories=33,saveat=bg.x_grid,reltol=reltol,dense=false) 
@time gsol = solve(eprob,Tsit5(),EnsembleGPUArray(),trajectories=33,saveat=bg.x_grid,reltol=reltol,dense=false) 
println("solved gpu ensemble")

#serial solve
@time for i in length(k_grid)
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k_grid[i], ℓᵧ, ℓ_ν, ℓ_mν,n_q)
    prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (xᵢ , zero(T)), hierarchy)
    solve(prob,KenCarp4(),saveat=bg.x_grid,reltol=reltol,dense=false)(bg.x_grid)
end

