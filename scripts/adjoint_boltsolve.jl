# Basic attempt to use sciml
using OrdinaryDiffEq, SciMLSensitivity
using SimpleChains
using Random
rng = Xoshiro(123);
using Plots
# Bolt boilerplate
using Bolt
kMpch=0.01;
ℓᵧ=3;
reltol=1e-5;
abstol=1e-5;
p_dm=[log(0.3),log(- -3.0)]
Ω_c,α_c = exp(p_dm[1]),-exp(p_dm[2]) #log pos params
#standard code
𝕡 = CosmoParams{Float32}(Ω_c=Ω_c,α_c=α_c)
𝕡

bg = Background(𝕡; x_grid=-20.0f0:0.1f0:0.0f0)

typeof(Bolt.η(-20.f0,𝕡,zeros(Float32,5),zeros(Float32,5)))


ih= Bolt.get_saha_ih(𝕡, bg);
k = 𝕡.h*kMpch  #get k in our units

typeof(8π)

# hierarchytestnn = Hierarchy_nn(BasicNewtonian(), 𝕡test, bgtest, ihtest, kk[1], p1,15);

# hierarchytestspl = Hierarchy_spl(BasicNewtonian(), 𝕡test, bgtest, ihtest, kk[1], 
#                                  Bolt.spline(ones(length(bgtest.x_grid)),bgtest.x_grid),
#                                  Bolt.spline(ones(length(bgtest.x_grid)),bgtest.x_grid),
#                                  15);

# u0test = Bolt.initial_conditions_nn(-20.0,hierarchytestnn);

# hierarchy_nn!(zeros(pertlen+2), u0test, hierarchytestnn, -20.0 )

hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ);
results = boltsolve(hierarchy; reltol=reltol, abstol=abstol);
res = hcat(results.(bg.x_grid)...)
# δ_c,v_c = res[end-3,:],res[end-2,:]


# Deconstruct the boltsolve...
#-------------------------------------------
xᵢ = first(Float32.(hierarchy.bg.x_grid))
u₀ = Float32.(initial_conditions(xᵢ, hierarchy));
# prob = ODEProblem{true}(hierarchy_nn!, u₀, (xᵢ , zero(T)), hierarchy)
# sol = solve(prob, ode_alg, reltol=reltol, abstol=abstol,
#             saveat=hierarchy.bg.x_grid,  
#             )
u₀

# NN = get_nn(4,pertlen+2,32)
# nnin = hcat([u...,k,x])

# NN setup
m=16;
pertlen=2(ℓᵧ+1)+5
NN₁ = SimpleChain(static(pertlen+2),
    TurboDense{true}(tanh, m), 
    TurboDense{true}(tanh, m), 
    TurboDense{false}(identity, 2) #have not tested non-scalar output
    );
p1 = SimpleChains.init_params(NN₁;rng); 
G1 = SimpleChains.alloc_threaded_grad(NN₁); 



# plot the solution to see if jagged
Plots.plot(bg.x_grid,sol(bg.x_grid)[end-3,:])



Plots.plot!(bg.x_grid,results(bg.x_grid)[end-3,:])
Plots.plot!(bg.x_grid,solt(bg.x_grid)[end-3,:],ls=:dash)
Plots.plot(bg.x_grid,sol(bg.x_grid)[end-4,:])
Plots.plot!(bg.x_grid,results(bg.x_grid)[end-4,:])
Plots.plot!(bg.x_grid,solt(bg.x_grid)[end-4,:],ls=:dash)



function hierarchy_nn_p!(du, u, p, x; hierarchy=hierarchy,NN=NN₁)
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)
    α_c = par.α_c
    Θ, Θᵖ, Φ, δ_c, v_c,δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    Θ′, Θᵖ′, _, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6

    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]
                                #   + Ω_c * a^(4+α_c) * σ_c 
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ_c
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        )
    # matter
    nnin = hcat([u...,k,x])
    println(size(nnin))
    u′ = NN(nnin,p)
    println(size(u′))
    δ′ = u′[1] #k / ℋₓ * v - 3Φ′
    v′ = u′[2]

    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[1] + v_b)
    # photons
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
    #END RSA

    du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end


function hierarchy_nn_p(u, p, x; hierarchy=hierarchy,NN=NN₁)
    # compute cosmological quantities at time x, and do some unpacking
    k, ℓᵧ, par, bg, ih = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih
    Ω_r, Ω_b, Ω_c, H₀² = par.Ω_r, par.Ω_b, par.Ω_c, bg.H₀^2 
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)
    csb² = ih.csb²(x)
    α_c = par.α_c

    Θ = [u[1],u[2],u[3],u[4]]
    Θᵖ = [u[5],u[6],u[7],u[8]]
    Φ, δ_c, v_c,δ_b, v_b = u[9:end]
    # Θ, Θᵖ, Φ, δ_c, v_c,δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    # Θ′, Θᵖ′, _, _, _, _, _ = unpack(du, hierarchy)  # will be sweetened by .. syntax in 1.6

    # metric perturbations (00 and ij FRW Einstein eqns)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[3]
                                #   + Ω_c * a^(4+α_c) * σ_c 
                                  )

    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_c * a^(2+α_c) * δ_c
        + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[1]
        )
    # matter
    nnin = hcat([u...,k,x])
    u′ = NN(nnin,p)
    # δ′, v′ = NN(nnin,p)
    δ′ = u′[1] #k / ℋₓ * v - 3Φ′
    v′ = u′[2]

    δ_b′ = k / ℋₓ * v_b - 3Φ′
    v_b′ = -v_b - k / ℋₓ * ( Ψ + csb² *  δ_b) + τₓ′ * R * (3Θ[2] + v_b)
    # photons
    Π = Θ[3] + Θᵖ[3] + Θᵖ[1]

    # Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    # Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    Θ′0 = -k / ℋₓ * Θ[2] - Φ′
    Θ′1 = k / (3ℋₓ) * Θ[1] - 2k / (3ℋₓ) * Θ[3] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[2] + v_b/3)
    Θ′2 = 2 * k / ((2*2+1) * ℋₓ) * Θ[3-1] -
            (2+1) * k / ((2*2+1) * ℋₓ) * Θ[3+1] + τₓ′ * (Θ[3] - Π * δ_kron(2, 2) / 10)
    # for ℓ in 2:(ℓᵧ-1 )
        # Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
        #     (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    # end
    # polarized photons
    # Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    Θᵖ′0 = -k / ℋₓ * Θᵖ[2] + τₓ′ * (Θᵖ[1] - Π / 2)
    # for ℓ in 1:(ℓᵧ-1)
    Θᵖ′1 = 1 * k / ((2*1+1) * ℋₓ) * Θᵖ[2-1] -
        (1+1) * k / ((2*1+1) * ℋₓ) * Θᵖ[2+1] + τₓ′ * (Θᵖ[2] - Π * δ_kron(1, 2) / 10)
    Θᵖ′2 = 2 * k / ((2*2+1) * ℋₓ) * Θᵖ[3-1] -
        (2+1) * k / ((2*2+1) * ℋₓ) * Θᵖ[3+1] + τₓ′ * (Θᵖ[3] - Π * δ_kron(2, 2) / 10)
    # end
    # for ℓ in 1:(ℓᵧ-1)
    #     Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
    #         (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    # end
    # photon boundary conditions: diffusion damping
    # Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[ℓᵧ]
    # Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - ( (ℓᵧ + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[ℓᵧ]
    Θ′3 = k / ℋₓ * Θ[4-1] - ( (3 + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θ[4]
    Θᵖ′3 = k / ℋₓ * Θᵖ[4-1] - ( (3 + 1) / (ℋₓ * ηₓ) - τₓ′ ) * Θᵖ[4]
    # du[2(ℓᵧ+1)+1:2(ℓᵧ+1)+5] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    # println("type return: ", typeof([Θ′0, Θ′1, Θ′2, Θ′3, Θᵖ′0, Θᵖ′1, Θᵖ′2, Θᵖ′3 , Φ′, δ′, v′, δ_b′, v_b′  ]) )
    return [Θ′0, Θ′1, Θ′2, Θ′3, Θᵖ′0, Θᵖ′1, Θᵖ′2, Θᵖ′3 , Φ′, δ′, v′, δ_b′, v_b′  ]
end
#-------------------------------------------
# test the input
hierarchy_nn_p(u₀,p1,-20.0f0)
typeof(u₀)

du_test_p = rand(pertlen+2);
hierarchy_nn_p!(du_test_p,u₀,p1,-20.0)

du_test_p

u₀
prob = ODEProblem{false}(hierarchy_nn_p, u₀, (xᵢ , zero(Float32)), p1)
sol = solve(prob, KenCarp4(), reltol=reltol, abstol=abstol,
            saveat=hierarchy.bg.x_grid,  
            )
# attempt to solve

#now try again with the sciml sensitivity

# Write some loss functions (sciml g functions) in terms of P(k) for *single mode*
# But first do some basic tests
dg_dscr(out, u, p, t, i) = (out.=-1.0.+u)
ts = -20.0:1.0:0.0

# This is dLdP * dPdu bc dg is really dgdu
PL = Bolt.get_Pk # Make some kind of function like this extracting from existing Pk function by feeding in u.
Nmodes = # do the calculation for some fiducial box size (motivated by a survey)
σₙ = P ./ Nmodes
dPdL = -2.f0 * (data-P) ./ σₙ.^2 # basic diagonal loss
dPdu = # the thing to extract, some weighted ratio of the matter density perturbations, so only 3 elements of u
dg_dscr_P(out,u,p,t,i) =(out.= -1.0)



g_cts(u,p,t) = (sum(u).^2) ./ 2.f0
dg_cts(out, u, p, t) = (out .= 2.f0*sum(u))

# This seems like it takes forever even just to produce an error? Is it just the first time or is it because 
# it has to do finite diff?
# maybe it is too many ts? Maybe the network just sucks so the solver takes a long time...

# discrete
@time res = adjoint_sensitivities(sol,#results,
                            KenCarp4(),sensealg=InterpolatingAdjoint(autojacvec=false,autodiff=false);t=ts,dgdu_discrete=dg_dscr,abstol=abstol,
                            reltol=reltol); # works (or at least spits out something)
# after the first one, this time returns:
#326.239041 seconds (3.36 G allocations: 102.113 GiB, 7.42% gc time) - for full bg ts
#140.153200 seconds (1.46 G allocations: 44.512 GiB, 7.59% gc time) - for ts with spacing dx=1.0

# continuous
@time res = adjoint_sensitivities(sol,#results,
                            KenCarp4(),sensealg=InterpolatingAdjoint(autojacvec=true);dgdu_continuous=dg_cts,g=g_cts,abstol=abstol,
                            reltol=reltol); # 
# after the first one, this time returns:
#111.447511 seconds (1.24 G allocations: 37.651 GiB, 7.26% gc time)

# Alright great so this won't work until we can make Bolt take an array as input, even abstract/component array
# I bet Zack actually tried this previously and this is what he found? Can ask him...

# Idk if we can do this with the code because of the interpolators...
# What does hierarchy actually hold that is the problem?



@time res = adjoint_sensitivities(sol,
                            KenCarp4(),sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true));
                            dgdu_continuous=dg_cts,g=g_cts,abstol=abstol,
                            reltol=reltol); # 


"""
it holds some ints for maximum boltzmann hierarchy, obviously this is not a problem for compnent arrays

Also what are we even talking about? We just want to make it recognize the NN parameters as p right?

So we should consider hierarchy etc. fixed, and it should be no problem to solve the backwards thing with that struct
We can just make it global?

Didn't I already encounter this in May-ish? What did I do there? I think I redefined boltsolve
Clearly the boltsolve here is NOT RIGHT

Ok but even the one I defined to extend Hierarchy to hold NN parameters won't work in this sense bc still holds interpolators

Can we just hack the code to make Hierarchy struct in to an "AbstractArray"? 
Ask Zack/Marius...wait for an answer

Alright, so, what I'll do for now is just make this file hold a global set of cosmological parameters.
Then we can freeze a global version of the background and ionization history interpolators
    Later, As Zack said, we can, without too much trouble, extend the code to just take a component array of Chebyshev weights,
    and then construct our interpolators out of those every time instead of using the splines.
    (Or we could even use the bsplines, we just save the bases and the weights, do what I did for the du/du' thing for every query of it)
    We just lose the overall functionality and ease of use we've already relied on Julia packages for...
After that those no longer have to be passed to the hierarchy! function, the only thing that will have to be passed is NN weights
At that point this SciML machinery should just work...
So I basically just need to redefine hierarchy! to only take the nn parameters and maek sure it can access the global interpolators
(In principle we could just run HMC over cosmo params this way by managing the scope a little better, but would be gross and probably slow)

This should not be hard -> 30 min exercise...

Ok this kind of works.

Now: 
1. Update the cost functions to be something we actually care about (power spectrum in a single k mode, generalize to multiple later)
2. Try the autodiff vecs (which will probably break) -> it does
3. Rough performance numbers (profile) what is taking so long? -> performannce btwn 2 is similar, reducing dscrt pts by factor of 10 is only performance gain of 2
4. Embed this in an optimization framework (make it a black box gradient function call) and test that gradient
5. Run that optimization framework (which may have to be on nersc, which is fine)
6. Is this actually faster than doing forwarddiff for the gradient? It should be...in principle
7. Can we *make* it faster than forwarddiff gradient? (Tolerance relaxation over training, fewer discrete points etc.)

"""

# Does pretraining with "good" parameters make this faster, since it is smooth?

#check what the derivative function actually looks like for this solution
#^It looks right early on when dm doesn't matter so much, then later on in MD looks super wrong
#^Basically what we expect

# Try pre-training on \LambdaCDM
# 1. Generate f(u,t) from true solution u
# 2. Train NN on f(u,t) (w/ e.g. SC train_unbatched)
# 3. Now evaluate the derivative with adjoints and compare the times to before
# (if it is much faster this is an argument to enforce smoothness somehow)
#^or maybe we just initialize the weights to be numerically small, but not zero?
aa = zeros(Float64,pertlen)
hierarchy_nn_p!(aa,results(ts[1]),p1,ts[1])
aa
# train SC network for many steps
fu_test = zeros(Float32,length(bg.x_grid),pertlen);
for i in 1:length(bg.x_grid)
    fu_test[i,:] .= hierarchy_nn_p(results(bg.x_grid[i]),p1,bg.x_grid[i])
    # println(fu_test[i,:])
end

fu_test_dc_vc = hcat(fu_test[:,end-3],fu_test[:,end-2]);
fu_test_dc_vc


scloss = SimpleChains.add_loss(NN₁, SquaredLoss(fu_test_dc_vc))
λ=1e-1
scloss_reg = FrontLastPenalty(scloss,
                                    L2Penalty(λ), L2Penalty(λ) ) 
N_ADAM = 100_000
typeof(results(ts))
# X_train = hcat(Array(results(ts))', ts, k.*ones(length(ts)))
X_train = hcat(Array(results(bg.x_grid))', k.*ones(length(bg.x_grid)),bg.x_grid)
size(X_train')

NN₁(X_train[1,:],p1)
NN₁(X_train',p1)
scloss(X_train[1,:],p1)
scloss(X_train',p1)
scloss(X_train,p1)
X_train'
SimpleChains.train_unbatched!(
                G1, p1, scloss,
                X_train', SimpleChains.ADAM(), N_ADAM 
            );
p1
scloss(X_train',p1)



solt = solve(remake(prob), KenCarp4(), u0 = u₀, tspan = (-20.0, 0.0), p =  p1 
            );
