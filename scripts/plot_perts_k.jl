using Revise
using Bolt
using ForwardDiff
#using PyPlot
using Plots
pyplot()
using BenchmarkTools
using Printf
#using ThreadPools

#input ingredients
par = CosmoParams()
bg = Background(par;x_grid=-20.0:0.1:0.0)#,logq_grid=-8.0:0.1:1.0)
ih = IonizationHistory(Peebles(), par, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100) #quadratically spaced k points
#numerical parameters

#look at the f0 and dlnf0, they look fine
plot(bg.logq_grid,10.0 .^(2. *bg.logq_grid) .* bg.f0)
plot!(collect(-6:0.2:1),(10.0 .^ collect(-6:0.2:1)) .^2 .* bg.f0(collect(-6:0.2:1)),ls=:dash) #check interpolant looks okay, crudely yes
plot(bg.logq_grid,log10.(abs.(bg.df0)))

ℓᵧ=8 #cutoff
ℓ_ν=10 #not used except for size here, should pass
ℓ_mν=ℓ_ν
n_q=10
reltol=1e-5 #cheap  rtol
#solve hierarchy at single x to check
x=-8 #just picking a number
a=exp(x)
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
println("pert vector length=",pertlen)
results=zeros(pertlen,length(k_grid))
#do single k for simplic
hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k_grid[1], ℓᵧ)
perturb = boltsolve(hierarchy; reltol=reltol)
u = perturb(x)
println(u)

#@btime @qthreads
#seems to take like a minute (without threads) (rel nus)
#with mass nus now takes around (~4-5mins) at much reduced reltol (1e-11-1e-5) with dumb integration
for (i_k, k) in enumerate(k_grid)
    println(i_k)
    hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ)
    perturb = boltsolve(hierarchy; reltol=reltol)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_k] = u #z should use unpack somehow
end

#Integrate the q moments for ℳ0 and ℳ2 for plotting
ℳρ,ℳσ = zeros(length(k_grid)),zeros(length(k_grid))
#bgℳρ,bgℳσ = ρP_0(a,par) #actually don't think we need this
for (i_k, k) in enumerate(k_grid)
    ℳρ[i_k],ℳσ[i_k] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_k],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_k],
                            bg,exp(x),par) .*(a^(-4) /bg.ρ_crit )
end

#'Need to figure out the units of k here'
labels = [raw"$\Phi$",raw"$\delta$",raw"$v$",raw"$|\delta_{b}|$",raw"$|v_{b}|$"]
plot(legend=:bottomleft)
for i in 1:5
    plot!(log10.(k_grid/ bg.H₀),
          log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+i,:])),
          label=labels[i])
end
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[1,:])),label=raw"$|\Theta_{0}|$")

plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+1,:])),
      label=raw"$|\mathcal{N}_{0}|$")
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+3,:])),
      label=raw"$|\mathcal{N}_{2}|$")

plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+1,:])),
      label=raw"$|\mathcal{M}_{0q0}|$")
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+2,:])),
    label=raw"$|\mathcal{M}_{0q2}|$")

#these look okay except for normalization...
plot!(log10.(k_grid/ bg.H₀),log10.(abs.(ℳρ)),
      label=raw"$|\mathcal{M}_{0}|$")
plot!(log10.(k_grid/ bg.H₀),log10.(abs.(ℳσ)),
      label=raw"$|\mathcal{M}_{2}|$")

plot(log10.(k_grid/ bg.H₀),1 ./(ℳρ./(results[2(ℓᵧ+1)+1,:])))
plot(log10.(k_grid/ bg.H₀),1 ./(ℳσ./(results[2(ℓᵧ+1)+3,:])))
ylims!(0,1e3)
# plot ℓ=1 perts:
# plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2,:])),label=raw"$|\Theta_{1}|$")
# plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+2,:])),
#     label=raw"$|\mathcal{N}_{1}|$")
#xlims!(1,3)

ylabel!(raw"$\delta_{i}(k)$")
xlabel!(raw"$k \ H₀$")
title!("z=$(@sprintf("%.0f", exp(-x)-1))")
savefig("../compare/bolt_perts_k_z$(@sprintf("%.0f", exp(-x)-1)).png")
