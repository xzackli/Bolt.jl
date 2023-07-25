# Try replacing Lux with SimpleChains
seed = parse(Int64, ARGS[1])
println("seed = $(seed)")

import Pkg
Pkg.activate("/global/cfs/cdirs/m4051/ude/Bolt.jl")
using OrdinaryDiffEq
using Random
using LinearAlgebra, Statistics,Lux
using Bolt
using Plots
using Optimization,SciMLSensitivity,OptimizationOptimisers,ComponentArrays,OptimizationOptimJL
using AbstractDifferentiation
import AbstractDifferentiation as AD, ForwardDiff
using DelimitedFiles
rng = Xoshiro(seed)


println("Done loading packages")
# setup ks won't use all of them here...
L=2f3
lkmi,lkmax,nk = log10(2.0f0*π/L),log10(0.2f0),8
kk = 10.0f0.^(collect(lkmi:(lkmax-lkmi)/(nk-1):lkmax))
ℓᵧ=15
pertlen=2(ℓᵧ+1)+5

# define network
m = 12
U = Lux.Chain(Lux.Dense(pertlen+1, m, tanh), #input is u,t
                  Lux.Dense(m, m,tanh),
                  Lux.Dense(m, 2))
p, st = Lux.setup(rng, U)



# copy the hierarchy function os it works for nn - for this to work you need the Hierarchy_nn struct and unpack_nn in perturbations.jl
function hierarchy_nn!(du, u, hierarchy::Hierarchy_nn{T, BasicNewtonian}, x)  where T
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)

    # the new free cdm index (not used here)
    α_c = par.α_c
    # get the nn params
    p_nn = hierarchy.p

    Θ, Θᵖ, Φ, δ_c, v_c,δ_b, v_b = unpack_nn(u, hierarchy)  
    Θ′, Θᵖ′, _, _, _, _, _ = unpack_nn(du, hierarchy) 

    # Here I am throwing away the neutrinos entriely, which is probably bad
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                #   + Ω_c * a^(4+α_c) * σ_c  #ignore this
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ_c
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        )

    # matter
    nnin = hcat([u...,x])
    û = U(nnin,p_nn,st)[1]
    δ′ = û[1]
    v′ = û[2]
    # here we implicitly assume σ_c = 0

    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)
    # photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    for ℓ in 2:(ℓᵧ-1)
        Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * Bolt.δ_kron(ℓ, 2) / 10)
    end
    Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓᵧ-1)
        Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[ℓᵧ]
    du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′
    return nothing
end

# use only the longest k mode
function hierarchy_nnu!(du, u, p, x) 
    hierarchy = Hierarchy_nn(BasicNewtonian(), 𝕡test, bgtest, ihtest, kk[1],p,15);
    hierarchy_nn!(du, u, hierarchy, x)
end



# some setup
tspan = (-20.0f0, 0.0f0)
𝕡test = CosmoParams(Ω_c=0.3,α_c=-3.0);
bgtest = Background(𝕡test; x_grid=-20.0f0:1f-1:0.0f0);
ihtest = Bolt.get_saha_ih(𝕡test, bgtest);
hierarchytest = Hierarchy(BasicNewtonian(), 𝕡test, bgtest, ihtest, kk[1],15);
hierarchytestnn = Hierarchy_nn(BasicNewtonian(), 𝕡test, bgtest, ihtest, kk[1], ComponentArray(p),15);
u0 = Bolt.initial_conditions_nn(tspan[1],hierarchytestnn)

# problem for truth and one we will remake
prob_trueode = ODEProblem(Bolt.hierarchy!, u0, tspan,hierarchytest)
prob_nn = ODEProblem(hierarchy_nnu!, u0, (bgtest.x_grid[1],bgtest.x_grid[end]), ComponentArray(p))

# Generate some noisy data (at reduced solver tolerance)
ode_data = Array(solve(prob_trueode, KenCarp4(), saveat = bgtest.x_grid,
                abstol = 1e-3, reltol = 1e-3))
δ_true,v_true = ode_data[end-3,:],ode_data[end-2,:]
σfakeode = 0.1
noise_fakeode = δ_true .*  randn(rng,size(δ_true)).*σfakeode
noise_fakeode_v = v_true .*  randn(rng,size(v_true)).*σfakeode
Ytrain_ode = hcat([δ_true .+ noise_fakeode,v_true .+ noise_fakeode_v]...)
# noise_both = Float32.(hcat([noise_fakeode,noise_fakeode_v]...))
noise_both = Float32.(hcat([δ_true.*σfakeode,v_true.*σfakeode]...))

# NB I dropped the "Float64" type argument to the ComponentArray, maybe we should put it back?
function predict(θ, T = bgtest.x_grid)
    _prob = remake(prob_nn, u0 = u0, tspan = (T[1], T[end]), p =  θ)
    res = Array(solve(_prob, KenCarp4(), saveat = T,
                abstol = 1e-3, reltol = 1e-3))
    return hcat(res[end-3,:],res[end-2,:])
end

#log loss
# function loss(θ)
#     X̂ = predict(θ)
#     log(mean(abs2, (Ytrain_ode - X̂)./noise_both ) )
# end
#raw loss
function loss(θ)
    X̂ = predict(θ)
    log(mean(abs2, (Ytrain_ode - X̂)) )
end

# adtype = Optimization.AutoZygote()
adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, ComponentVector(p))


# Training
losses,pp = [],[];
callback = function (p, l)
    push!(losses, l)
    push!(pp, p)
    if length(losses) % 5 == 0
    println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    return false
end
niter_0,η_0 = 20,1.0
niter_1,niter_2,niter_3 =50,50,20
η_1,η_2,η_3 = 1.0,0.1,0.01
res1 = Optimization.solve(optprob, ADAM(η_1), callback = callback, maxiters = niter_1)
test_predict_o1 = predict(res1.u)

optprob2 = remake(optprob,u0 = res1.u);
res2 = Optimization.solve(optprob2, ADAM(η_2), callback = callback, maxiters = niter_2)
test_predict_o2 = predict(res2.u)

optprob3 = remake(optprob2,u0 = res2.u);
res3 = Optimization.solve(optprob3, ADAM(η_3), callback = callback, maxiters = niter_3)
test_predict_o3 = predict(res3.u)

println("done optimizing")

# save the results
writedlm("../../data/t_raw_enn_v1_seed$(seed).dat",predict(res3.u))
writedlm("../../data/t_raw_loss_v1_seed$(seed).dat",losses)
writedlm("../../data/t_raw_params_v1_seed$(seed).dat",pp)

println("done saving")

function get_Φ′_Ψ(u,hierarchy::Hierarchy{T},x) where T
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)
    α_c = par.α_c
    Θ, Θᵖ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  
    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ 
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        )

    return Φ′,Ψ
end

# The reconstructed function
Φ′_true,Ψ_true = zeros(length(hierarchytest.bg.x_grid)),zeros(length(hierarchytest.bg.x_grid))
nn_δ′,nn_v′ = zeros(length(hierarchytest.bg.x_grid)),zeros(length(hierarchytest.bg.x_grid))
for j in 1:length(hierarchytest.bg.x_grid)
    Φ′_true[j],Ψ_true[j] = get_Φ′_Ψ(ode_data[:,j],hierarchytest,hierarchytest.bg.x_grid[j])
    nnin = hcat([ode_data[:,j]...,hierarchytest.bg.x_grid[j]])
    nn_u′ = U(nnin,res3.u,st)[1]
    nn_δ′[j],nn_v′[j] = nn_u′[1], nn_u′[2]
end


true_δ′ = hierarchytestnn.k ./ hierarchytestnn.bg.ℋ .* v_true .- 3Φ′_true
true_v′ = -v_true .- hierarchytestnn.k  ./ hierarchytestnn.bg.ℋ  .* Ψ_true

writedlm("../../data/t_raw_true_v1_seed$(seed).dat",[true_δ′,true_v′])
writedlm("../../data/t_raw_recon_v1_seed$(seed).dat",[nn_δ′,nn_v′])