using Revise
using Bolt
using ForwardDiff
#using PyPlot
using Plots
pyplot()
using BenchmarkTools
using Printf
using QuadGK
#using ThreadPools

#input ingredients
par = CosmoParams()
logqmin,logqmax = -6,-1
n_q = 15
logq_pts = logqmin:(logqmax-logqmin)/(n_q-1):logqmax
bg = Background(par;x_grid=-20.0:0.1:0.0,nq=n_q) #,logq_grid=logq_pts)
ih = IonizationHistory(Peebles(), par, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100) #quadratically spaced k points

# #quickly check bg against nr approx
# Tγ = (15/ π^2 *bg.ρ_crit *par.Ω_r)^(1/4)
# νfac = (90 * 1.2020569 /(11 * π^4)) * (par.Ω_r * par.h^2 / Tγ)#the factor that goes into nr approx to neutrino energy density
# par.Σm_ν*νfac/par.h^2 *((par.N_ν/3)^(3/4)) /(ρP_0(1,par)[1]/bg.ρ_crit)
# (par.N_ν/3)^(3/4)

function f00(q)
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 *3.9669896e-11 *5.042e-5)^(1/4) ##assume instant decouple for now
    gs =  2 #should be 2 for EACH neutrino family (mass eigenstate)
    return gs / (2π)^3 / ( exp(q/Tν) +1)
end

function dlnf0dlnq0(q) #this is actually only used in perts
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 *3.9669896e-11 *5.042e-5)^(1/4) ##assume instant decouple for now
    return -q / Tν /(1 + exp(-q/Tν))
end

#find correct factor for normalization...
ρν0 = 7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r * bg.ρ_crit / 2 #used to be div by Neff, now div by floor(Neff-1)
ρν=ρν0*(exp(-20))^(-4)
aaa=4π  * quadgk(q ->  q^2 * -dlnf0dlnq0(q) *q * f00(q),
            1e-6, 1e-1,rtol=1e-6)[1]/4/ρν0
#!
aaa
# Removed the splines from bg so these won't work anymore
# #check the splining error:
# #use both splines - error of ~1.8e-3
# 4π  * quadgk(q ->  q^2 * -bg.df0(log10(q)) *q * bg.f0(log10(q)),
#             1e-6, 1e-1,rtol=1e-6)[1]/4/ρν0
#
# #only use df0 spline - error of ~ 5e-6, ~ rtol
# 4π  * quadgk(q ->  q^2 * -bg.df0(log10(q)) *q * f00(q),
#             1e-6, 1e-1,rtol=1e-6)[1]/4/ρν0
#
# #only use f0 spline - again we get error of ~1.8e-3 - so this is the problem...
# 4π  * quadgk(q ->  q^2 * -dlnf0dlnq0(q) *q * bg.f0(log10(q)),
#             1e-6, 1e-1,rtol=1e-6)[1]/4/ρν0

#test that ρ_σ is the same as bg when passed ones - it is up to quadgk tol...
bgrho,_ =  (exp(-20)^(-4)) .* ρ_σ(ones(n_q) ,
               zeros(n_q),bg,exp(-20),par)
ρP_0(exp(-20),par)
ρν #analytic answer

#@btime @qthreads
ℓᵧ=8 #cutoff
ℓ_ν=10 #not used except for size here, should pass
ℓ_mν=ℓ_ν
reltol=1e-5 #cheaper  rtol
#solve hierarchy at single x to check
x=-8 #just picking a number
a=exp(x)
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
println("pert vector length=",pertlen)
results=zeros(pertlen,length(k_grid))
n_q
for (i_k, k) in enumerate(k_grid)
    println(i_k)
    hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ, n_q)
    perturb = boltsolve(hierarchy; reltol=reltol)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_k] = u #z should use unpack somehow
end
results

#Integrate the q moments for ℳ0 and ℳ2 for plotting
ℳρ,ℳσ = zeros(length(k_grid)),zeros(length(k_grid))
ℳθ = zeros(length(k_grid))
Ω_ν =  7*(2/3)*par.N_ν/8 *(4/11)^(4/3) *par.Ω_r
norm𝒩 = 1/(4Ω_ν * bg.ρ_crit / 2)
for (i_k, k) in enumerate(k_grid)
    ℳρ[i_k],ℳσ[i_k] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_k],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_k],
                            bg,a,par)#.*norm𝒩
    ℳθ[i_k],_ = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q,i_k],
                            zeros(n_q),
                            bg,a,par)

end
labels = [raw"$\Phi$",raw"$\delta$",raw"$v$",raw"$|\delta_{b}|$",raw"$|v_{b}|$"]
plot(legend=:bottomleft)
#title!("ICs")
for i in 1:5
    plot!(log10.(k_grid/ bg.H₀),
          log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+i,:])),
          label=labels[i])
end
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[1,:])),label=raw"$|\Theta_{0}|$")

plot!(log10.(k_grid/ bg.H₀), log10.(abs.( results[2(ℓᵧ+1)+1,:])),
      label=raw"$|\mathcal{N}_{0}|$")
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+2,:])),
      label=raw"$|\mathcal{N}_{1}|$")
plot!(log10.(k_grid/ bg.H₀), log10.(abs.(results[2(ℓᵧ+1)+3,:])),
      label=raw"$|\mathcal{N}_{2}|$")

results[2(ℓᵧ+1)+1,:] ./ ℳρ
results[2(ℓᵧ+1)+2,:] ./ ℳθ
results[2(ℓᵧ+1)+3,:] ./ ℳσ

plot!(log10.(k_grid/ bg.H₀),log10.(abs.(ℳρ)),
      label=raw"$|\mathcal{M}_{0}|$",ls=:dash)
plot!(log10.(k_grid/ bg.H₀),log10.(abs.(ℳθ)),
    label=raw"$|\mathcal{M}_{1}|$",ls=:dash)
plot!(log10.(k_grid/ bg.H₀),log10.(abs.(ℳσ)),
      label=raw"$|\mathcal{M}_{2}|$",ls=:dash)


ylabel!(raw"$\delta_{i}(k)$")
xlabel!(raw"$k / H0$")
title!("z=$(@sprintf("%.0f", exp(-x)-1))")
savefig("../compare/bolt_perts_k_z$(@sprintf("%.0f", exp(-x)-1)).png")
