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
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
x_grid = collect(-20:0.01:0.0)

# Choose a k-mode to compare to saved class perturbations at
k_options = ["p03", "p3", "1p0", #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
             "p01", ] #newly computed k modes
k_choice = k_options[1]
#Read in CLASS perturbations
#CLASS keys (for reference):
#['k (h/Mpc)', 'd_g', 'd_b', 'd_cdm', 'd_ur', 'd_ncdm[0]', 'd_tot',
#'phi', 'psi', 't_g', 't_b', 't_cdm', 't_ur', 't_ncdm[0]', 't_tot']
retnf = open( @sprintf("./test/data/zack_N_class_px_k%s_nofluid_nonu.dat",k_choice),"r" ) do datafile
# an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#the second column is just a repeated k value, so remember it and delete col
kclass = retnf[2][1] #read class k mode from file (in h/Mpc)
# k = (bg.H₀*3e5/100)*kclass #get k in our units ->old value
k = (bg.H₀*299792.458/100)*kclass #get k in our units
class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))

xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))] #horizon crossing ish
println("k = ", kclass," log10k = ", log10(kclass), " h/Mpc")

#pert setup
ℓᵧ=50
ℓ_ν=50
ℓ_mν=20
reltol=1e-6
abstol=1e-6
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid)) #arrays for the massive neutrino integrated perts
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
#solve (with rsa)
results_with_rsa = boltsolve_rsa(hierarchy; reltol=reltol, abstol=abstol)
# results_with_rsa = boltsolve(hierarchy; reltol=reltol, abstol=abstol)


# Residual plot against CLASS
#reverse the arrays because ow complains
ϵ=1e-2
x_grid[1:end-1].-class_pxsnf[1,:][end:-1:1][1]
min_class_idx = minimum(findall(<(ϵ), x_grid[1:end-1]./class_pxsnf[1,:][end:-1:1][1] .- (1-ϵ) ))
x_grid_aligned = x_grid[1:end-1][min_class_idx:end]
ours_2_Mpcm1 = 1/(2.1331e-35 *3e5) #unit conversion 1/([km/s/Mpc]*[c/km/s])


#checking for early super-horizon evoluion from ICs internal to bolt...
#the fact that this doesn't use rsa doesn't matter here
# plot(x_grid,results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]/results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,1],
#      )
# plot!(x_grid,results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]/results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,1],
#      )
# plot!(x_grid,results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]/results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,1],
#      )
# plot!(x_grid,results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]/results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,1],
#      )
# hline!([1],ls=:dot,color=:black)
# xlims!(-20,-13)
# ylims!(0.999,1.001)
# xlabel!(raw"$x$")
# ylabel!(raw"$\delta_{i}(x)/\delta_{i}(x=-20)$")


##Compare ratios of perturbation evolution with CLASS

#phi
itpphiclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[7,:][end:-1:1])
plot(x_grid_aligned, (results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:])[1:end-1][min_class_idx:end]./itpphiclass.(x_grid_aligned), label="phi" )
hline!([1],ls=:dot,color=:black,label=false )

#matter density
itpdelclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[4,:][end:-1:1])
plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:])[1:end-1][min_class_idx:end]./itpdelclass.(x_grid_aligned), label="mat" )

#baryon density
itpbarclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[3,:][end:-1:1])
plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:])[1:end-1][min_class_idx:end]./itpbarclass.(x_grid_aligned), label="bar")
# plot(x_grid_aligned, itpbarclass.(x_grid_aligned)./itpbarclass_rf.(x_grid_aligned), label="bar class hyrec/rf" )

#photon monopole
# itpgamclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[2,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end]./itpgamclass.(x_grid_aligned), label="pho" )

#massless neutrino monopole
# itpnu0class = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[5,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+1,:]* 𝕡.h*4)[1:end-1][min_class_idx:end]./itpnu0class.(x_grid_aligned), label="nu0" )
# vline!([this_rsa_switch],ls=:dot,color=:green,label="rsa")
# vline!([log(1/1101)],ls=:dot,color=:green,label="~decoupling")
# vline!([log(1/3400)],ls=:dot,color=:gray,label="~MD")


vline!([xhor], color="red")
ylims!(1.0 - 0.2e-2, 1.0 + 0.2e-2)
# xlims!(-8,0)
xlabel!(raw"$x$")
ylabel!(raw"$\frac{\delta_{i,Bolt}}{\delta_{i,CLASS}}(x)$")
title!("Compare CLASS - Bolt - k=$(@sprintf("%.3f", kclass))")
# savefig("../compare/614nonu_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

##
plot(x_grid_aligned, (results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end])
plot!(x_grid_aligned, -itpgamclass.(x_grid_aligned), label="CLASS", ls=:dash)
plot!(ylim=(1.5,2.0))

##

# plot(x_grid_aligned, itpphiclass.(x_grid_aligned), label="pho" )
# plot(x_grid_aligned, -(results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end] .- itpgamclass.(x_grid_aligned), label="pho" )

##

# itpbarclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[3,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:])[1:end-1][min_class_idx:end]./itpbarclass.(x_grid_aligned), label="bar")

# plot(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:])[1:end-1][min_class_idx:end])
# plot!(x_grid_aligned, itpbarclass.(x_grid_aligned), label="CLASS" )

