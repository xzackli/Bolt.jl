using Revise
using Bolt
using Plots
using Printf
using Interpolations

𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
logqmin,logqmax = -6,-1
x_grid = collect(-10:0.1:0)

ℓᵧ=1000
ℓ_ν=100
ℓ_mν=50
reltol=1e-5 #cheaper  rtol
k =  1000bg.H₀*.3/.333 /10
kbolt = k/(bg.H₀*3e5/100)
xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))] #horizon crossing

println("k = ", kbolt,
        " log10k = ", log10(kbolt), " h/Mpc")
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5

results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid))

hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
perturb = boltsolve(hierarchy; reltol=reltol)

for (i_x, x) in enumerate(x_grid)
    println(i_x)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
    ℳρ[i_x],ℳσ[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
                            bg,exp(x),𝕡)
    #normalization for plotting, divide by integral of just momentum measure
    # ℳρ[i_x]=ℳρ[i_x] ./ bg.ρ₀ℳ(x)
    ℳρ[i_x]=ℳρ[i_x] ./ (ρ_σ(ones(length(bg.quad_pts)),
                                   zeros(length(bg.quad_pts)),
                                   bg,exp(x),𝕡)[1] )

    println(ρP_0(exp(x),𝕡,bg.quad_pts,bg.quad_pts)[1],' ', bg.ρ₀ℳ(x),' ',(exp(x)^-4 *ρ_σ(ones(length(bg.quad_pts)),
                                   zeros(length(bg.quad_pts)),
                                   bg,exp(x),𝕡)[1] ))
end
results

#What?? Why should ρP_0 be different from ρ₀ℳ which is computed using the same function??
ρP_0(1.0,𝕡,bg.quad_pts,bg.quad_pts)[1]
ρ_σ(ones(length(bg.quad_pts)),
                               zeros(length(bg.quad_pts)),
                               bg,1,𝕡)[1]

bg.ρ₀ℳ(0)

#CLASS perturbations
#CLASS keys:
#['k (h/Mpc)', 'd_g', 'd_b', 'd_cdm', 'd_ur', 'd_ncdm[0]', 'd_tot',
#'phi', 'psi', 't_g', 't_b', 't_cdm', 't_ur', 't_ncdm[0]', 't_tot']
# ret = open("./test/data/class_px_kp3.dat","r") do datafile
ret = open("./test/data/class_px_kp03.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#By default CLASS uses fluid approximation, which introduces almost 2x error for massive neutrinos at lower x
#We don't want to compare to this
# retnf = open("./test/data/class_px_kp3_nofluid.dat","r") do datafile
retnf = open("./test/data/class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end

#the second column is just a repeated k value, so remember it and delete col
kclass = ret[2][1]
class_pxs = transpose(reduce(hcat,ret[1:end .!= 2]))
class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))
class_pxs
println("kclass is ", kclass, " kbolt is ",kbolt, " ratio (c/b) is ", kclass/kbolt)

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
plot(class_pxs[1,:],log10.(abs.(class_pxs[8,:])),
    label=raw"$\Phi_{\rm{CLASS}}$")
plot!(x_grid, log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h),
      label=raw"$h \Phi_{\rm{Bolt}}$",ls=:dash)
vline!([xhor],ls=:dot,c=:black,label="k/[2πℋ(x)h]=1")


xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt (NR) - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/nr_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

#massless neutrino monopole 𝒩0
plot(class_pxs[1,:],log10.(abs.(class_pxs[5,:])),
     label=raw"$\nu_{0,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid, log10.(abs.(4results[2(ℓᵧ+1)+1,:]* 𝕡.h)),
      label=raw"$4h \nu_{0,\rm{Bolt}}$",ls=:dash)

#photon Θ0 monopole
plot(class_pxs[1,:],log10.(abs.(class_pxs[2,:])),
      label=raw"$\Theta_{0,\rm{CLASS}}$",legend=:topleft)
plot!(x_grid, log10.(abs.(results[1,:]* 𝕡.h*4)),
      label=raw"$4 h \Theta_{0,\rm{Bolt}}$",ls=:dash)

#massive neutrino monopole ℳ0
plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[6,:])),
    label=raw"$m\nu_{0,\rm{CLASS,nf}}$",
    legend=:topleft)
plot!(class_pxs[1,:],log10.(abs.(class_pxs[6,:])),
    label=raw"$m\nu_{0,\rm{CLASS,f}}$",
    ls=:dot)
plot!(x_grid, log10.(abs.(ℳρ* 𝕡.h)),
    label=raw"$h m\nu_{0,\rm{Bolt}}$",ls=:dash)
    #ls=:dot)
vline!([xhor],ls=:dot,c=:black,label=raw"$k/(2\pi a H h)=1$")

xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt (R) - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/r_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

#Look at some ratios with CLASS
#reverse the arrays because ow complains
itpnuclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[6,:][end:-1:1])
itpnuclassf = LinearInterpolation(class_pxs[1,:][end:-1:1],class_pxs[6,:][end:-1:1])
#drop the last bolt element because arrays are strangely aligned...
plot(x_grid[1:end-1], ((-ℳρ* 𝕡.h)[1:end-1]./itpnuclass.(x_grid[1:end-1]) ))
plot!(x_grid[1:end-1], ((-ℳρ* 𝕡.h)[1:end-1]./itpnuclassf.(x_grid[1:end-1]) ))
println(typeof(class_pxs[1,:]), ' ', typeof(class_pxs[6,:]))
hline!([1],ls=:dot,color=:black)
#fluid approx looks similar in shape to 1104.2935

#check Phi, delta
itpphiclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[8,:][end:-1:1])
plot!(x_grid[1:end-1], (results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h )[1:end-1]./itpphiclass.(x_grid[1:end-1]) )
hline!([1],ls=:dot,color=:black)

itpdelclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[4,:][end:-1:1])
plot(x_grid[1:end-1], -(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h )[1:end-1]./itpdelclass.(x_grid[1:end-1]) )
hline!([1],ls=:dot,color=:black)

itpbarclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[3,:][end:-1:1])
plot(x_grid[1:end-1], -(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h )[1:end-1]./itpbarclass.(x_grid[1:end-1]) )
hline!([1],ls=:dot,color=:black)
