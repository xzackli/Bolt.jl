using OrdinaryDiffEq
using Random
using LinearAlgebra, Statistics,Lux
rng = Xoshiro(123)
using Bolt
using Plots
using Optimization,SciMLSensitivity,OptimizationOptimisers,ComponentArrays,OptimizationOptimJL
using AbstractDifferentiation
import AbstractDifferentiation as AD, ForwardDiff
using SimpleChains
using AdvancedHMC,LogDensityProblems

# setup ks won't use all of them here...
L=2f3
lkmi,lkmax,nk = log10(2.0f0*π/L),log10(0.2f0),8
kk = 10.0f0.^(collect(lkmi:(lkmax-lkmi)/(nk-1):lkmax))
ℓᵧ=15
pertlen=2(ℓᵧ+1)+5

# define network
width=8
# U = Lux.Chain(Lux.Dense(pertlen+1, width, tanh), #input is u,t
#                   Lux.Dense(width, width,tanh),
#                   Lux.Dense(width, 2))
# p, st = Lux.setup(rng, U)
function get_nn(m,d_in)
    NN = SimpleChain(static(d_in),
        TurboDense{true}(tanh, m), 
        TurboDense{true}(tanh, m), 
        TurboDense{false}(identity, 2) #have not tested non-scalar output
    );
    p = SimpleChains.init_params(NN;rng); 
    G = SimpleChains.alloc_threaded_grad(NN);
    return NN,p,G
end

NN,p,G = get_nn(width,pertlen+1)


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
    # û = U(nnin,p_nn,st)[1]
    # û = sc_noloss(nnin,p_nn)'  .* std(true_u′,dims=1) .+ mean(true_u′,dims=1)
    û = NN(nnin,p_nn)' 
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
hierarchytestnn = Hierarchy_nn(BasicNewtonian(), 𝕡test, bgtest, ihtest, kk[1],p,15);
u0 = Bolt.initial_conditions_nn(tspan[1],hierarchytestnn)

NN(hcat([u0...,-20.f0]),p)
# NN([-20.f0],p)

# problem for truth and one we will remake
prob_trueode = ODEProblem(Bolt.hierarchy!, u0, tspan,hierarchytest)
# prob_nn = ODEProblem(hierarchy_nnu!, u0, (bgtest.x_grid[1],bgtest.x_grid[end]), ComponentArray{Float64}(p))

# Generate some noisy data (at reduced solver tolerance)
ode_data = Array(solve(prob_trueode, KenCarp4(), saveat = bgtest.x_grid,
                abstol = 1f-3, reltol = 1f-3))
δ_true,v_true = ode_data[end-3,:],ode_data[end-2,:]
σfakeode = 0.1f0
noise_fakeode = δ_true .*  randn(rng,size(δ_true)).*σfakeode
noise_fakeode_v = v_true .*  randn(rng,size(v_true)).*σfakeode
Ytrain_ode = hcat([δ_true .+ noise_fakeode,v_true .+ noise_fakeode_v]...)
noise_both_old = Float32.(hcat([noise_fakeode,noise_fakeode_v]...))
noise_both = Float32.(hcat([δ_true.*σfakeode,v_true.*σfakeode]...))

#float conversion
fl_xgrid = Float32.(bgtest.x_grid)
fu0 = Float32.(u0)
prob_nn = ODEProblem(hierarchy_nnu!, fu0, (fl_xgrid[1],fl_xgrid[end]), p)

function predict(θ, T = fl_xgrid)
    _prob = remake(prob_nn, u0 = fu0, tspan = (T[1], T[end]), p =  θ)
    res = Array(solve(_prob, KenCarp4(), saveat = T,
                abstol = 1f-3, reltol = 1f-3))
    return hcat(res[end-3,:],res[end-2,:])
end

#log loss
function loss(θ)
    X̂ = predict(θ)
    log(sum(abs2, (Ytrain_ode - X̂)./ noise_both ) )
end


# adtype = Optimization.AutoZygote()
adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf,p)



# Training
function make_gif_plot(prediction,iter)
    p = Plots.plot()
    Plots.scatter!(p,bgtest.x_grid,Ytrain_ode[:,1],label="data")
    Plots.plot!(p,bgtest.x_grid,prediction,lw=2,label="iter = $(iter)")
    Plots.xlabel!(p,raw"$\log(a)$")
    Plots.ylabel!(p,raw"$\delta_{c}(a)$")
    savefig("../plots/learning_v1_multnoise$(σfakeode)_Adam$(iter)_s123_sc_chisq_to.png")
    p
end

losses = [];
pp =[];
callback = function (p, l)
    push!(losses, l)
    push!(pp, p)
    if length(losses) % 5 == 0
    println("Current loss after $(length(losses)) iterations: $(losses[end])")
    make_gif_plot(predict(p)[:,1],length(losses))
    end
    return false
end
niter_1,niter_2,niter_3 =50,50,20
η_1,η_2,η_3 = 1.f0,0.1f0,0.01f0
#loss doesn't go down very much at all with 1f-3 η1

res1 = Optimization.solve(optprob, ADAM(0.01), callback = callback, maxiters = niter_1)
# this is pretty slow, on the order of half an hour, but it is at least running!
#Current loss after 50 iterations: 3.8709838260084415
test_predict_o1 = predict(res1.u)

optprob2 = remake(optprob,u0 = res1.u);
res2 = Optimization.solve(optprob2, ADAM(η_2), callback = callback, maxiters = niter_2)
test_predict_o2 = predict(res2.u)
#Current loss after 100 iterations: 1.1319215191941459

optprob3 = remake(optprob,u0 = res2.u);
res3 = Optimization.solve(optprob3, ADAM(η_3), callback = callback, maxiters = 100)
test_predict_o3 = predict(res3.u)
#Current loss after 120 iterations: 1.0862281116475199

# FIXME MORE INVOLVED ADAM SCHEDULE
# for the heck of it try BFGS, not sure what parameters it usually takes
optprob4 = remake(optprob3,u0 = res4.u);
# res4 = Optimization.solve(optprob4, BFGS(), 
res4 = Optimization.solve(optprob4, ADAM(η_3/1000), callback = callback, 
                            maxiters = 100)
                            # callback = callback, maxiters = 10)
test_predict_o4 = predict(res4.u)
# wow this is slow...I guess due to Hessian approximation?
# We have the gradient so idk why this takes so much longer than ADAM?
# Somehow loss actually goes up also? Maybe we overshoot, can try again with specified initial stepsize...
test_predict_o1[:,1]
# Plots of the learned perturbations
Plots.scatter(bgtest.x_grid,Ytrain_ode[:,1],label="data")
Plots.plot!(bgtest.x_grid,δ_true,label="truth",lw=2.5)#,yscale=:log10)
Plots.plot!(bgtest.x_grid,predict(pc)[:,1],label="opt-v1",lw=2.5,ls=:dash)
Plots.title!(raw"$\delta_{c}$")
Plots.xlabel!(raw"$\log(a)$")
Plots.ylabel!(raw"$\delta_{c}(a)$")
savefig("../plots/deltac_learning_v1_multnoise$(σfakeode)_Adam$(niter_1)_$(niter_2)_$(niter_3)_$(η_1)_$(η_2)_$(η_3)_bfgs.png")
pc
p
log10.(test_predict_o4[:,2])


# It *SEEMS LIKE* the model isn't flexible enough - i.e. when I shift to
# weighted loss from square loss the late exponential part gets worse
# while the early part gets worse.
# For SC training earlier, regularization helped a bit but not much...

Plots.scatter(bgtest.x_grid,Ytrain_ode[:,2],label="data")#,legend=:bottomright)
Plots.plot!(bgtest.x_grid,v_true,label="truth")
Plots.plot!(bgtest.x_grid,test_predict_o1[:,2],label="opt-v1")
Plots.plot!(bgtest.x_grid,test_predict_o4[:,2],label="opt-v1-full",lw=2.5,ls=:dot)
Plots.plot!(bgtest.x_grid,predict(pc)[:,2],label="opt-v1",lw=2.5,ls=:dash)
Plots.title!(raw"$v_{c}$")
Plots.xlabel!(raw"$\log(a)$")
Plots.ylabel!(raw"$v_{c}(a)$")
savefig("../plots/vc_learning_v1_multnoise$(σfakeode)_Adam$(niter_1)_$(niter_2)_$(niter_3)_$(η_1)_$(η_2)_$(η_3)_bfgs.png")

