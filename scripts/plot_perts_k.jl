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
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 1000) #quadratically spaced k points
println(bg.H₀) #the units on this are in H0?
 #must need to account for evoln somehow?
#numerical parameters
ℓᵧ=8 #cutoff
ℓ_ν=10 #not used except for size here, should pass
reltol=1e-11

#solve hierarchy at single x to check
x=-8 #just picking a number
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+5
println("pert vector length=",pertlen)
results=zeros(pertlen,length(k_grid))

#@btime @qthreads
#seems to take like a minute (without threads)
for (i_k, k) in enumerate(k_grid)
    hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ)
    perturb = boltsolve(hierarchy; reltol=reltol)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_k] = u #z should use unpack somehow
end
#'Need to figure out the units of k here'
labels = [raw"$\Phi$",raw"$\delta$",raw"$v$",raw"$|\delta_{b}|$",raw"$|v_{b}|$"]
plot()
for i in 1:5
    plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+i,:])),label=labels[i])
end
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[1,:])),label=raw"$|\Theta_{0}|$")

plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+1,:])),
      label=raw"$|\mathcal{N}_{0}|$")
# plot ℓ=1 perts:
# plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2,:])),label=raw"$|\Theta_{1}|$")
# plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+2,:])),
#     label=raw"$|\mathcal{N}_{1}|$")
#xlims!(1,3)

ylabel!(raw"$\delta_{i}(k)$")
xlabel!(raw"$k \ H₀$")
title!("z=$(@sprintf("%.0f", exp(-x)-1))")
savefig("../compare/bolt_perts_k_z$(@sprintf("%.0f", exp(-x)-1)).png")
