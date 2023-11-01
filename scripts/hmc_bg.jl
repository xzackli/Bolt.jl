using SimpleChains
using Plots,CairoMakie,PairPlots
using LinearAlgebra
using AdvancedHMC, ForwardDiff
using OrdinaryDiffEq
using LogDensityProblems
using BenchmarkTools
using DelimitedFiles
using Plots.PlotMeasures
using DataInterpolations
using Random,Distributions
rng = Xoshiro(123)

# function to get a simplechain (fixed network architecture)
function train_initial_network(Xtrain,Ytrain,rng; λ=0.0f0, m=128,use_l_bias=true,use_L_bias=false,N_ADAM=1_000,N_epochs = 3,N_rounds=4)
    #regularized network
    d_in,d_out = size(Xtrain)[1],size(Ytrain)[1]
    mlpd = SimpleChain(
        static(d_in),
        TurboDense{use_l_bias}(tanh, m), 
        TurboDense{use_l_bias}(tanh, m), 
        TurboDense{use_L_bias}(identity, d_out) #have not tested non-scalar output
        )
    p = SimpleChains.init_params(mlpd;rng); 
    G = SimpleChains.alloc_threaded_grad(mlpd); 
    mlpdloss = SimpleChains.add_loss(mlpd, SquaredLoss(Ytrain))
    mlpdloss_reg = FrontLastPenalty(SimpleChains.add_loss(mlpd, SquaredLoss(Ytrain)),
                                    L2Penalty(λ), L2Penalty(λ) ) 
    loss = λ > 0.0f0 ?  mlpdloss_reg : mlpdloss
    for k in 1:N_rounds
        for _ in 1:N_epochs #FIXME, if not looking at this, don't need to split it up like so...
            SimpleChains.train_unbatched!(
                G, p, loss, Xtrain, SimpleChains.ADAM(), N_ADAM 
            );
        end
    end
    mlpd_noloss = SimpleChains.remove_loss(mlpd)
    return mlpd_noloss,p
end


function run_hmc(ℓπ,initial_θ;n_samples=2_000,n_adapts=1_000,backend=ForwardDiff)
    D=size(initial_θ)
    metric = DiagEuclideanMetric(D)#ones(Float64,D).*0.01)
    hamiltonian = Hamiltonian(metric, ℓπ, backend)
    initial_ϵ = find_good_stepsize(hamiltonian, initial_θ)
    integrator = Leapfrog(initial_ϵ)
    proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    samples, stats = sample(hamiltonian, proposal, initial_θ, n_samples, adaptor, n_adapts; progress=true);
    h_samples = hcat(samples...)'
    return h_samples = hcat(samples...)',stats
end

using Bolt
# first try just sampling with ODE over the 2 variables of interest with HMC and true ODE
function dm_model(p_dm::Array{DT,1};kMpch=0.01,ℓᵧ=15,reltol=1e-5,abstol=1e-5) where DT #k in h/Mpc
    Ω_c,α_c = exp(p_dm[1]),-exp(p_dm[2]) #log pos params
    # println("Ω_c,α_c = $(Ω_c) $(α_c)")
    #standard code
    𝕡 = CosmoParams{DT}(Ω_c=Ω_c,α_c=α_c)
    bg = Background(𝕡; x_grid=-20.0:0.1:0.0)
    # 𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    # ih = IonizationHistory(𝕣, 𝕡, bg)
    ih= Bolt.get_saha_ih(𝕡, bg);

    k = 𝕡.h*kMpch  #get k in our units
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ)
    results = boltsolve(hierarchy; reltol=reltol, abstol=abstol)
    res = hcat(results.(bg.x_grid)...)
    δ_c,v_c = res[end-3,:],res[end-2,:]
    return δ_c#,v_c
end

dm_model([log(0.3),log(- -3.0)])

# δ_true,v_true = dm_model([log(0.3),log(- -3.0)]);
δ_true= dm_model([log(0.3),log(- -3.0)]);
δ_true

#test jacobian
jdmtest = ForwardDiff.jacobian(dm_model,[log(0.3),log(- -3.0)])
jdmtest

#HMC for 2 ode params 
# this is obviously very artificial, multiplicative noise
σfakeode = 0.1
noise_fakeode = δ_true .*  randn(rng,size(δ_true)).*σfakeode
# noise_fakeode_v = v_true .*  randn(rng,size(v_true)).*σfakeode
# Ytrain_ode = [δ_true .+ noise_fakeode,v_true .+ noise_fakeode_v]
Ytrain_ode = δ_true .+ noise_fakeode
noise_fakeode
Plots.plot(δ_true)
Plots.scatter!(Ytrain_ode)

# density functions for HMC
struct LogTargetDensity_ode
    dim::Int
end

LogDensityProblems.logdensity(p::LogTargetDensity_ode, θ) = -sum(abs2, (Ytrain_ode-dm_model(θ))./noise_fakeode) / 2  # standard multivariate normal
LogDensityProblems.dimension(p::LogTargetDensity_ode) = p.dim
LogDensityProblems.capabilities(::Type{LogTargetDensity_ode}) = LogDensityProblems.LogDensityOrder{0}()
ℓπode = LogTargetDensity_ode(2)

initial_ode = [log(0.3),log(- -3.0)];
ode_samples, ode_stats = run_hmc(ℓπode,initial_ode;n_samples=200,n_adapts=100)

ode_samples

ode_stats
ode_labels = Dict(
 :α => "parameter 1",
 :β => "parameter 2"
)

# Annoyingly PairPlots provides very nice functionality only if you use DataFrames (or Tables)
α,β = ode_samples[:,1], ode_samples[:,2]
using DataFrames
df = DataFrame(;α,β)

pairplot(df  ,PairPlots.Truth( 
            (;α =log(0.3),β=log(- -3.0)), 
            label="Truth" )  )

LogDensityProblems.logdensity(ℓπode,[log(0.3),log(- -5.0)])

ForwardDiff.jacobian()


@btime LogDensityProblems.logdensity(ℓπode,initial_ode)
#  14.466 ms (7036 allocations: 1.73 MiB)
@profview LogDensityProblems.logdensity(ℓπode,initial_ode)
# so this is just slow, maybe not surprising

# open("./test/ode_deltac_hmc_trueinit_samplesshort.dat", "w") do io
#     writedlm(io, ode_samples)
# end
