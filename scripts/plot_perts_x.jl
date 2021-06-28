using Revise
using Bolt
using Plots
using Printf
using Interpolations
using DelimitedFiles

𝕡 = CosmoParams()
n_q=15
logqmin,logqmax = -6,-1
bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
x_grid = collect(-20:0.01:0.0)
# SOUND SPEED TESTING
# plot(x_grid,log10.(ih.csb²(x_grid)))
# writedlm("../compare/x_csb2_test.df", (x_grid,ih.csb²(x_grid)))


#check RSA condition on τ′
ih.g̃(-11)
plot(x_grid,ih.g̃(x_grid))

#Why a factor of conformal H is necessary? need to read abou this...
plot(x_grid,log10.(-ih.τ′.(x_grid) .* bg.η.(x_grid) .* bg.H₀))# .* bg.ℋ.(x_grid)))
plot!(x_grid,log10.(-ih.τ′.(x_grid) .* bg.η.(x_grid) .* bg.ℋ.(x_grid)))
plot!(x_grid,log10.(k .* bg.η.(x_grid) ),ls=:dash)
hline!([log10(5)])
hline!([log10(45)],ls=:dash)

#---
#see above plot but for this particular mode at this cosmology the k condition happens later
this_rsa_switch = x_grid[argmin(abs.(k .* bg.η.(x_grid) .- 45))]
#^This seems kinda late?

ℓᵧ=50
ℓ_ν=50
ℓ_mν=10
reltol=1e-5 #cheaper  rtol
k =  1000bg.H₀*.3/.333 /10
#fudge factors to match the CLASS k mode values
fudge_class_03 = 0.029158189805725376/0.030030030030030026 #fudge factor so we get the right .03
# fudge_class_3 = 0.3000962432008842/0.3003003003003003 #fudge factor so we get the right .3
fudge_class_lowres_03 = 0.020881445483105634/0.029158189805725376
kbolt = k/(bg.H₀*3e5/100)*fudge_class_03#*fudge_class_lowres_03

xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))] #horizon crossing ish

println("k = ", kbolt,
        " log10k = ", log10(kbolt), " h/Mpc")
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5

results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid))

hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
perturb = boltsolve(hierarchy; reltol=reltol)

# plot(log10.(abs.(perturb.u[1])))
perturb.destats

for (i_x, x) in enumerate(x_grid)
    # println(i_x,', ',x)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
    ℳρ[i_x],ℳσ[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
                            bg,exp(x),𝕡)
    #normalization for plotting, divide by integral of just momentum measure
    ℳρ[i_x]=ℳρ[i_x] ./ bg.ρ₀ℳ(x) .* exp(-4x) #missing factor of a^-4 from pert rho
    # ℳρ[i_x]=ℳρ[i_x] ./ (ρ_σ(ones(length(bg.quad_pts)),
    #                                zeros(length(bg.quad_pts)),
    #                                bg,exp(x),𝕡)[1] )

    # println(ρP_0(exp(x),𝕡,bg.quad_pts,bg.quad_pts)[1],' ', bg.ρ₀ℳ(x),' ',(exp(x)^-4 *ρ_σ(ones(length(bg.quad_pts)),
    #                                zeros(length(bg.quad_pts)),
    #                                bg,exp(x),𝕡)[1] ))
end

#THESE SHOULD BE THE SAME IN RSA, BUT DIFFEQ.jl NOT COOPERATING
results[2(ℓᵧ+1)+1,end]
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,end]
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
# goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end

#the second column is just a repeated k value, so remember it and delete col
kclass = retnf[2][1]
class_pxs = transpose(reduce(hcat,ret[1:end .!= 2]))
class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))
class_pxs
println("kclass is ", kclass, " kbolt is ",kbolt, " ratio (c/b) is ", kclass/kbolt)


#quick look at these - copying similar syntax from plot perts k
#skipping velocities this time just for simplicity
#matter δ
class_pxsnf[1,:],class_pxs[1,:]
plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[4,:])),
     label=raw"$\delta_{c,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid,log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h),
      label=raw"$h \delta_{\rm{Bolt}}$",ls=:dash)

#baryon δ_b
plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[3,:])),
    label=raw"$\delta_{b,\rm{CLASS}}$")
plot!(x_grid,log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h)),
    label=raw"$h \delta_{b,\rm{Bolt}}$",ls=:dash)
xlims!(-12,0)

#pick a time (x=-15) and check if Phi is equal or not
class_pxsnf[1,end-8]
x_grid[500]
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,1]* 𝕡.h
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,500]* 𝕡.h
class_pxsnf[8,end-8]
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,1], class_pxsnf[8,end-8] ./ 𝕡.h
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,1]* 𝕡.h ./class_pxsnf[8,end-8]
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,500], class_pxsnf[8,end-8] ./ 𝕡.h
results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,500]* 𝕡.h ./class_pxsnf[8,end-8]

#Check the initial conditions gaianst what CLASS does
#everything here is in units of k (not h/Mpc)
xini=-20
aini=exp(xini) #for now
#at initial time ignore neutrin mss, assume we got temp right
rhor = 𝕡.Ω_r*(1+(3/3)*(7𝕡.N_ν/8)*(4/11)^(4/3)) * aini^(-4)
rhom = (𝕡.Ω_m + 𝕡.Ω_b ) * aini^(-3)
denom = 1 + rhom/rhor
fν = 𝕡.Ω_r*(3/3)*(7𝕡.N_ν/8)*(4/11)^(4/3)  * aini^(-4) / rhor
fγ = 𝕡.Ω_r*1*aini^(-4) / rhor
fc = 𝕡.Ω_m *  aini^(-3) / rhom
fb = 𝕡.Ω_b *  aini^(-3) / rhom
δγ,δν = 4results[1,1],4results[2(ℓᵧ+1)+1,1]#multiply by 3/2?
ℳρ[1]
δb,δc = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,1],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,1]
δTot = ( fγ*δγ + fν*δν  +(denom-1)*( fb*δb + fc*δc ) )/denom
θγ,θν = results[2,1] * 3k,results[2(ℓᵧ+1)+2,1] *3k #really we are converting to velocity
θb = - results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5,1] #really this is velocity
# vTot = -( 4/3*( fγ*θγ + fν*θν ) +(denom-1)*fb*θb )/denom #k??
# 3bg.ℋ(xini)/(k) *vTot
# (bg.ℋ(xini)/k)^2
# # vTot = ( ( fγ*results[2,1] + fν*results[2(ℓᵧ+1)+2,1] ) -(denom-1)*fb/4*results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5,1] )/denom #k??
# ΦCL = -3/2 * (bg.ℋ(xini)/k)^2 * ( δTot + 3bg.ℋ(xini) * vTot )
# ΦCL/4
#^Something is wrong with velocities here


#for Φinit=1.0
fν
Ψinit = -1.0/(1+2/5 * fν)
𝕡.Ω_r
7*(3/3)*𝕡.N_ν/8 *(4/11)^(4/3) *𝕡.Ω_r
7*(3/3)*𝕡.N_ν/8 *(4/11)^(4/3)

#Plot Phi and Psi from class
plot(class_pxsnf[1,:],class_pxsnf[8,:]/𝕡.h,
    label=raw"$\Phi_{\rm{CLASS}}$")
plot!(class_pxsnf[1,:],class_pxsnf[9,:]/𝕡.h,
        label=raw"$\Psi_{\rm{CLASS}}$")
hline!([1],ls=:dash,label=false,c="black")
hline!([0.8594],ls=:dot,label=false,c="blue")

#throw in space metric Φ also
plot(class_pxsnf[1,:],log10.(class_pxsnf[8,:]),
    label=raw"$\Phi_{\rm{CLASS}}$")
plot!(x_grid, log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h  ),
      label=raw"$h \Phi_{\rm{Bolt}}$",ls=:dash)
vline!([xhor],ls=:dot,c=:black,label="k/[2πℋ(x)h]=1")
title!("ad hoc factor of .95 Phi")
savefig("../compare/adhoc_superhorizonp95factor.png")

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

vline!([this_rsa_switch],label="RSA switch",ls=:dot)
plot!(x_grid, log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h*4),
      label=raw"$4h \Phi_{\rm{Bolt}}$",ls=:dash)

#massless neutrino dipole 𝒩1
plot(x_grid, log10.(abs.(4results[2(ℓᵧ+1)+2,:]* 𝕡.h)),
      label=raw"$4h \nu_{0,\rm{Bolt}}$",ls=:dash)
vline!([this_rsa_switch],label="RSA switch",ls=:dot)
# doesn't work, don't have phi' output plot!(x_grid, log10.(abs.(-2*results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h*4)),
#       label=raw"$|h *-2\Phi'_{\rm{Bolt}}/k|$",ls=:dash)


#photon Θ0 monopole
plot(class_pxs[1,:],log10.(abs.(class_pxs[2,:])),
      label=raw"$\Theta_{0,\rm{CLASS}}$",legend=:topleft)
plot!(x_grid, log10.(abs.(results[1,:]* 𝕡.h*4)),
      label=raw"$4 h \Theta_{0,\rm{Bolt}}$",ls=:dash)
vline!([this_rsa_switch],label="RSA switch",ls=:dot)
plot!(x_grid, log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:] .+
                     results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5,:] ./ k .*ih.τ′.(x_grid).*bg.ℋ.(x_grid))*𝕡.h*4),
      label=raw"$4h (\Phi_{\rm{Bolt}}+ 1/k τₓ′ * v_b)$",ls=:dash)
ylims!(-.2,1)

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

#fluid approx looks similar in shape to 1104.2935


#check Phi, delta
itpphiclass = LinearInterpolation(class_pxsnf[1,:][end-8:-1:1],class_pxsnf[8,:][end-8:-1:1])
# plot(x_grid[500:end-1],itpphiclass.(x_grid[500:end-1]))
plot(class_pxsnf[1,:],class_pxsnf[8,:],label=raw"$\Phi_{\rm{CLASS}}$")
plot!(x_grid, results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]*𝕡.h,label=raw"$\Phi_{\rm{Bolt}}$",ls=:dash)
hline!([𝕡.h],ls=:dash,color=:black,label=false )
# ylims!(.695,.705)
vline!([x_grid[75]],ls=:dot,color=:black,label="step 75")
title!("new ics all species,mnu=0.06,neff=3.046 include massive 00")
xlabel!(raw"$x$")
ylabel!(raw"$\Phi$")
vline!([xhor],ls=:dot,c=:blue,label="k/[2πℋ(x)h]=1")
savefig("../compare/newics_includemassive_bgpt_phiproblem_mnu0.06_neff3.046.png")
#matter density
itpdelclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[4,:][end:-1:1])
plot(x_grid[1:end-1], -(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h )[1:end-1]./itpdelclass.(x_grid[1:end-1]), label="mat" )
hline!([1],ls=:dot,color=:black,label=false )

#baryon density
itpbarclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[3,:][end:-1:1])
plot!(x_grid[1:end-1], -(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h )[1:end-1]./itpbarclass.(x_grid[1:end-1]), label="bar")
hline!([1],ls=:dot,color=:black,label=false )

#photon monopole
itpgamclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[2,:][end:-1:1])
plot!(x_grid[1:end-1], -(results[1,:]* 𝕡.h*4)[1:end-1]./itpgamclass.(x_grid[1:end-1]), label="pho" )
hline!([1],ls=:dot,color=:black,label=false )

#massless neutrino monopole
itpnu0class = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[5,:][end:-1:1])
plot!(x_grid[1:end-1], -(results[2(ℓᵧ+1)+1,:]* 𝕡.h*4)[1:end-1]./itpnu0class.(x_grid[1:end-1]), label="nu0" )
hline!([1],ls=:dot,color=:black,label=false )
