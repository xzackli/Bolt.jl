# using Revise
using Bolt
using Plots
using Printf
using Interpolations
using DelimitedFiles

# bg/ion setup
𝕡 = CosmoParams()
n_q=15
logqmin,logqmax = -6,-1
bg = Background(𝕡; x_grid=-20.0:0.004:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
x_grid = collect(-20:0.004:0.0)

# Choose a k-mode to compare to saved class perturbations at
k_options = ["p03", "p3", "1p0"] #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
k_choice = k_options[1]
#Read in CLASS perturbations
#CLASS keys (for reference):
#['k (h/Mpc)', 'd_g', 'd_b', 'd_cdm', 'd_ur', 'd_ncdm[0]', 'd_tot',
#'phi', 'psi', 't_g', 't_b', 't_cdm', 't_ur', 't_ncdm[0]', 't_tot']
#reading files for two diff k modes
ret = open( @sprintf("./test/data/class_px_k%s.dat",k_choice),"r" ) do datafile #note that these should only be used for comparing nonreion
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#By default CLASS uses fluid approximation, which introduces almost 2x error for massive neutrinos at lower x
#So compare to no fluid case to see if hierarchy is right
retnf = open( @sprintf("./test/data/class_px_k%s_nofluid_re.dat",k_choice),"r" ) do datafile
# an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#the second column is just a repeated k value, so remember it and delete col
kclass = retnf[2][1] #read class k mode from file (in h/Mpc)
k = (bg.H₀*3e5/100)*kclass #get k in our units
class_pxs = transpose(reduce(hcat,ret[1:end .!= 2]))
class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))
dipole_fac = 2kclass #for later normalization

#see above plot but for this particular mode at this cosmology the k condition happens later
this_rsa_switch = x_grid[argmin(abs.(k .* bg.η.(x_grid) .- 45))]


xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))] #horizon crossing ish
println("k = ", kclass," log10k = ", log10(kclass), " h/Mpc")

#pert setup
ℓᵧ=50
ℓ_ν=50
ℓ_mν=20
reltol=1e-8 #cheaper  rtol
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid)) #arrays for the massive neutrino integrated perts
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
#solve
perturb = boltsolve(hierarchy; reltol=reltol)

#get results and compute massive neutrino integrated moments
for (i_x, x) in enumerate(x_grid)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
    ℳρ[i_x],ℳσ[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
                            bg,exp(x),𝕡)
    #normalization for plotting, divide by integral of just momentum measure
    ℳρ[i_x]=ℳρ[i_x] ./ bg.ρ₀ℳ(x) .* exp(-4x) #missing factor of a^-4 from pert rho
end

#I don't know how to do the DE.jl splining over all the perts at once so this is just an array...
results_with_rsa = boltsolve_rsa(hierarchy; reltol=reltol)


#Look at the massless species to check RSA

#photon Θ0 monopole
plot(class_pxs[1,:],log10.(abs.(class_pxsnf[2,:])),
      label=raw"$\Theta_{0,\rm{CLASS}}$",legend=:topleft)
plot!(x_grid, log10.(abs.(results_with_rsa[1,:]*4)),
      label=raw"$4 \Theta_{0,\rm{Bolt,rsa}}$",ls=:dash)
plot!(x_grid, log10.(abs.(results[1,:]*4)),
      label=raw"$4 \Theta_{0,\rm{Bolt}}$",ls=:dash)
vline!([this_rsa_switch],label="RSA switch",ls=:dot)
xlims!(-8,0)
xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
#photon Θ1 dipole
plot(class_pxs[1,:],log10.(abs.(class_pxsnf[10,:])),
     label=raw"$\Theta_{1,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid, log10.(abs.(results_with_rsa[2,:] * dipole_fac)),
      label=raw"$4 \Theta_{1,\rm{Bolt,rsa}}$",ls=:dash)
plot!(x_grid, log10.(abs.(results[2,:] * dipole_fac)),
      label=raw"$4 \Theta_{1,\rm{Bolt}}$",ls=:dash)
vline!([this_rsa_switch],label="RSA switch",ls=:dot)
xlims!(-5,1)
xlims!(-7,0)

#neutrino 𝒩0 monopole
plot(class_pxs[1,:],log10.(abs.(class_pxsnf[5,:])),
      label=raw"$\mathcal{N}_{0,\rm{CLASS}}$",legend=:topleft)
plot!(x_grid, log10.(abs.(results_with_rsa[1+2(ℓᵧ+1),:]*4)),
      label=raw"$4 \mathcal{N}_{0,\rm{Bolt}}$",ls=:dash)
vline!([this_rsa_switch],label="RSA switch",ls=:dot)
ylims!(-.2,1)
xlims!(-7,0)
#neutrino 𝒩1 dipole
plot(class_pxs[1,:],log10.(abs.(class_pxs[13,:])),
      label=raw"$\mathcal{N}_{0,\rm{CLASS}}$",legend=:topleft)
plot!(x_grid, log10.(abs.(results_with_rsa[2(ℓᵧ+1)+2,:]* dipole_fac)),
      label=raw"$4 \mathcal{N}_{0,\rm{Bolt}}$",ls=:dash)
vline!([this_rsa_switch],label="RSA switch",ls=:dot)
xlims!(-5,1)


# Residual plot against CLASS
#reverse the arrays because ow complains
ϵ=1e-2
x_grid[1:end-1].-class_pxsnf[1,:][end:-1:1][1]
min_class_idx = minimum(findall(<(ϵ), x_grid[1:end-1]./class_pxsnf[1,:][end:-1:1][1] .- (1-ϵ) ))
x_grid_aligned = x_grid[1:end-1][min_class_idx:end]
ours_2_Mpcm1 = 1/(2.1331e-35 *3e5) #unit conversion 1/([km/s/Mpc]*[c/km/s])

# Check massive species for PL

#matter δ
plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[4,:])),
     label=raw"$\delta_{c,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid,log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]),
      label=raw"$h \delta_{\rm{Bolt}}$",ls=:dash)

#baryon δ_b
plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[3,:])),
    label=raw"$\delta_{b,\rm{CLASS}}$",legend=:topleft)
plot!(x_grid,log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:])),
    label=raw"$h \delta_{b,\rm{Bolt}}$",ls=:dash)
xlims!(-8,0)

#massive neutrino monopole ℳ0
plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[6,:])),
    label=raw"$m\nu_{0,\rm{CLASS,nf}}$",
    legend=:topleft)
plot!(class_pxs[1,:],log10.(abs.(class_pxs[6,:])),
    label=raw"$m\nu_{0,\rm{CLASS,f}}$",
    ls=:dot)
plot!(x_grid, log10.(abs.(ℳρ)),
    label=raw"$h m\nu_{0,\rm{Bolt}}$",ls=:dash)
    #ls=:dot)
vline!([xhor],ls=:dot,c=:black,label=raw"$k/(2\pi a H h)=1$")


#phi
itpphiclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[8,:][end:-1:1])
plot(x_grid_aligned, (results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]*𝕡.h)[1:end-1][min_class_idx:end]./itpphiclass.(x_grid_aligned), label="phi" )
hline!([1],ls=:dot,color=:black,label=false )

#matter density
itpdelclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[4,:][end:-1:1])
plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h )[1:end-1][min_class_idx:end]./itpdelclass.(x_grid_aligned), label="mat" )

#baryon density
itpbarclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[3,:][end:-1:1])
plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h )[1:end-1][min_class_idx:end]./itpbarclass.(x_grid_aligned), label="bar")
# plot(x_grid_aligned, itpbarclass.(x_grid_aligned)./itpbarclass_rf.(x_grid_aligned), label="bar class hyrec/rf" )

#photon monopole
itpgamclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[2,:][end:-1:1])
plot!(x_grid_aligned, -(results_with_rsa[1,:]* 𝕡.h*4)[1:end-1][min_class_idx:end]./itpgamclass.(x_grid_aligned), label="pho" )

#massless neutrino monopole
itpnu0class = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[5,:][end:-1:1])
plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+1,:]* 𝕡.h*4)[1:end-1][min_class_idx:end]./itpnu0class.(x_grid_aligned), label="nu0" )
vline!([this_rsa_switch],ls=:dot,color=:green,label="rsa")
vline!([log(1/1101)],ls=:dot,color=:red,label="decoupling")
ylims!(0.95,1.05)
xlims!(-8,0)
xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/reion_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")
