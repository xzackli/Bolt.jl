using Revise
using Bolt
using Plots
using Printf

𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡;x_grid=-20.0:0.1:0.0,nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
logqmin,logqmax = -6,-1
x_grid = collect(-10:0.1:0)

ℓᵧ=8
ℓ_ν=ℓᵧ
ℓ_mν=ℓ_ν
reltol=1e-5 #cheaper  rtol
k =  1000bg.H₀*.3/.333 /10
kbolt = k/(bg.H₀*3e5/100)
println("k = ", kbolt,
        " log10k = ", log10(kbolt), " h/Mpc")
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5

results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid))
for (i_x, x) in enumerate(x_grid)
    println(i_x)
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_mν,n_q)
    perturb = boltsolve(hierarchy; reltol=reltol)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
    ℳρ[i_x],ℳσ[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
                            bg,exp(x),𝕡)
end
results

#recompute this while testing
ℳρtest,ℳσtest = zeros(length(x_grid)),zeros(length(x_grid))
for (i_x, x) in enumerate(x_grid)
    ℳρtest[i_x],ℳσtest[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
                            bg,exp(-20),𝕡)
end

#CLASS perturbations
#CLASS keys:
#['k (h/Mpc)', 'd_g', 'd_b', 'd_cdm', 'd_ur', 'd_ncdm[0]', 'd_tot',
#'phi', 'psi', 't_g', 't_b', 't_cdm', 't_ur', 't_ncdm[0]', 't_tot']
# ret = open("./test/data/class_px_kp3.dat","r") do datafile
ret = open("./test/data/class_px_kp03.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end

#the second column is just a repeated k value, so remember it and delete col
kclass = ret[2][1]
class_pxs = transpose(reduce(hcat,ret[1:end .!= 2]))
class_pxs
println("kclass is ", kclass, " kbolt is ",kbolt)

#quick look at these - copying similar syntax from plot perts k
#skipping velocities this time just for simplicity
#matter δ
plot(class_pxs[1,:],log10.(abs.(class_pxs[4,:])),
     label=raw"$\delta_{c,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid,log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h),
      label=raw"$h \delta_{\rm{Bolt}}$",ls=:dash)

#baryon δ_b
plot!(class_pxs[1,:],log10.(abs.(class_pxs[3,:])),
    label=raw"$\delta_{b,\rm{CLASS}}$")
plot!(x_grid,log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h)),
    label=raw"$h \delta_{b,\rm{Bolt}}$",ls=:dash)

#throw in space metric Φ also
plot!(class_pxs[1,:],log10.(abs.(class_pxs[8,:])),
    label=raw"$\Phi_{\rm{CLASS}}$")
plot!(x_grid, log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h),
      label=raw"$h \Phi_{\rm{Bolt}}$",ls=:dash)


xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt (NR) - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/nr_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

#massless neutrino monopole 𝒩0
plot(class_pxs[1,:],log10.(abs.(class_pxs[5,:])),
     label=raw"$\nu_{0,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid, log10.(abs.(results[2(ℓᵧ+1)+1,:]* 𝕡.h*4)),
      label=raw"$4 h \nu_{0,\rm{Bolt}}$",ls=:dash)

#photon Θ0 monopole
plot!(class_pxs[1,:],log10.(abs.(class_pxs[2,:])),
      label=raw"$\Theta_{0,\rm{CLASS}}$")
plot!(x_grid, log10.(abs.(results[1,:]* 𝕡.h*4)),
      label=raw"$4 h \Theta_{0,\rm{Bolt}}$",ls=:dash)

#massive neutrino monopole ℳ0
plot!(class_pxs[1,:],log10.(abs.(class_pxs[6,:])),
    label=raw"$m\nu_{0,\rm{CLASS}}$",
    #legend=:topleft)
plot!(x_grid, log10.(abs.(ℳρ* 𝕡.h *4)),
    label=raw"$4 h m\nu_{0,\rm{Bolt}}$",ls=:dash)
plot!(x_grid, log10.(abs.(ℳρtest* 𝕡.h *4)),
    label=raw"$\mathrm{rel norm} 4 h m\nu_{0,\rm{Bolt}}$",ls=:dash)

xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt (mnu) - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/mnu_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")
