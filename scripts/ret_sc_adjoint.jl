using SimpleChains
using ForwardDiff
using OrdinaryDiffEq
using DataInterpolations
using Random,Distributions
rng = Xoshiro(123)
using Optimization,SciMLSensitivity,OptimizationOptimisers#,ComponentArrays # for fwd diff
using SciMLSensitivity

# Setup
#------------------------------------------------------------------------

# general NN function
function get_nn(m,d_in)
    NN = SimpleChain(static(d_in),
        TurboDense{true}(tanh, m), 
        # TurboDense{true}(tanh, m),  
        TurboDense{false}(identity, 1) 
    );
    p = SimpleChains.init_params(NN;rng); 
    G = SimpleChains.alloc_threaded_grad(NN);
    return NN,p,G
end

# Simple ODE problem
function f_true!(du,u,p,t) 
    du .= p.*u 
    return nothing 
end
u0 = [1.f0];
tspan = (0.f0,1.f0);
λp = -1.0f0;
prob_true = ODEProblem(f_true!,u0,tspan, λp)
tt = 1f-2:1f-2:1.f0
rtol=atol=1e-12

# Generate noisy "data" solution
σrel = 0.01f0; #multiplicative noise
true_soln = solve(prob_true,Tsit5(), saveat = tt,
                    abstol=atol,reltol=rtol)

#training data with noise
Ytrain = vcat(Array(true_soln)...).*(1.f0.+σrel.*randn(rng,eltype(u0),size(tt))) 
Ytrain_interp = CubicSpline(Ytrain,tt); #interpolated training data

# neural network function (here we don't let it depend on t bc it doesn't need it, but in general it should)
function f_nn!(du,u,p,t)
    du .=  NN(u,p)
    return nothing 
end

# NN ODE setup (initialization)
NN,p₀,G = get_nn(16,1); # shallow narrow network for testing
prob_nn = ODEProblem{true}(f_nn!,u0,tspan,p₀);
nn_soln = solve(prob_nn,Tsit5(),abstol=atol,reltol=rtol)

function predict_interp(θ)
    _prob = remake(prob_nn, u0 = u0, tspan = tspan, p =  θ)
        res = solve(_prob, Tsit5(), #saveat = collect(tt),
                abstol = atol, reltol = rtol,dense=true)
    return res
end

function loss_cts(θ;Yinterp=Ytrain_interp,ttfine = 1f-2:1f-3:1.f0) 
    X̂ = predict_interp(θ)(ttfine)
    sum((Yinterp(ttfine) .- Array(X̂)[1,:]).^2)*(ttfine[2]-ttfine[1]) #Riemann sum approx of cts integral
end

# ---> This is where I had a long training step, but took it out for simplicity (can put it back)

# Get the sciML adjoints
dg_cts(out, u, p, t) = (out .= -2.f0*(Ytrain_interp(t) - u[1]))
g_cts(u,p,t) = (Ytrain_interp(t) - u[1])^2
sciml_adj_soln_cts = adjoint_sensitivities(predict_interp(p₀),# (p is not learned, so  will be wrong)
                                    Tsit5(),
                                    g=g_cts, dgdu_continuous = dg_cts, 
                                    abstol = atol,reltol = rtol, # note in general we probably don't want to use the same fwd/back tols
                                    # iabstol=atol_sciml,ireltol=rtol_sciml
                                    sensealg = InterpolatingAdjoint(autojacvec=false, autodiff=true),#( autojacvec = true, autodiff=true)
                                    );

# This is the adjoint gradient (structure of the returne thing is (λ(t=t₀), grad) ):
sciml_adj_grad_cts = sciml_adj_soln_cts[2]


# Compare to FwdDiffGradient
fwd_diff_grad_cts = ForwardDiff.gradient(loss_cts,p₀)

prod( abs.(sciml_adj_grad_cts[1,:] .- fwd_diff_grad_cts) .<4f-5 )#true


# INSERT HMC CALLING THIS GRADIENT ON THE ABOVE LOSS HERE