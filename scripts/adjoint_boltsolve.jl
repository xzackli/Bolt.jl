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

hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ);
results = boltsolve(hierarchy; reltol=reltol, abstol=abstol);
res = hcat(results.(bg.x_grid)...)


# Deconstruct the boltsolve...
#-------------------------------------------
xᵢ = first(Float32.(hierarchy.bg.x_grid))
u₀ = Float32.(initial_conditions(xᵢ, hierarchy));

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

dg_dscr(out, u, p, t, i) = (out.=-1.0.+u)
ts = -20.0:1.0:0.0

# # This is dLdP * dPdu bc dg is really dgdu
# PL = Bolt.get_Pk # Make some kind of function like this extracting from existing Pk function by feeding in u.
# Nmodes = # do the calculation for some fiducial box size (motivated by a survey)
# σₙ = P ./ Nmodes
# dPdL = -2.f0 * (data-P) ./ σₙ.^2 # basic diagonal loss
# dPdu = # the thing to extract, some weighted ratio of the matter density perturbations, so only 3 elements of u
# dg_dscr_P(out,u,p,t,i) =(out.= -1.0)

g_cts(u,p,t) = (sum(u).^2) ./ 2.f0
dg_cts(out, u, p, t) = (out .= 2.f0*sum(u))

#

@time res = adjoint_sensitivities(sol,
                            KenCarp4(),sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true));
                            dgdu_continuous=dg_cts,g=g_cts,abstol=abstol,
                            reltol=reltol); # 
