# using Revise
using Bolt
using Plots
using Printf
using Interpolations, DataInterpolations
using DelimitedFiles

# bg/ion setup
𝕡 = CosmoParams()
n_q=15
logqmin,logqmax = -6,-1
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
ih = IonizationHistory(𝕣, 𝕡, bg)


# function customion(par, bg, Xₑ_function, Tmat_function, csb²_function)
#      x_grid = bg.x_grid
#      τ, τ′ = Bolt.τ_functions(x_grid, Xₑ_function, par, bg.ℋ)
#      g̃ = Bolt.g̃_function(τ, τ′)
#      spline, spline_∂ₓ, spline_∂ₓ² = Bolt.spline, Bolt.spline_∂ₓ, Bolt.spline_∂ₓ²
#      Xₑ_ = spline(Xₑ_function.(x_grid), x_grid)
#      τ_ = spline(τ.(x_grid), x_grid)
#      g̃_ = spline(g̃.(x_grid), x_grid)
#      Tmat_ = spline(Tmat_function.(x_grid), x_grid)
#      csb²_ = spline(csb²_function.(x_grid), x_grid)
 
#      return IonizationHistory(
#            (τ(0.)),
#          Xₑ_,
#          τ_,
#          spline_∂ₓ(τ_, x_grid),
#          spline_∂ₓ²(τ_, x_grid),
#          g̃_,
#          spline_∂ₓ(g̃_, x_grid),
#          spline_∂ₓ²(g̃_, x_grid),
#          Tmat_,
#            #spline_∂ₓ(Tmat_, x_grid),
#            csb²_,
#          # Trad_ #why do we even have this?
#      )
#  end
 
# retnf = open( @sprintf("./test/data/zack_N_class_px_nofluid_nonu_x_e.dat"),"r" ) do datafile
#      [parse.(Float64, split(line)) for line in eachline(datafile)]
# end

# class_a = retnf[1]
# class_x = log.(class_a)
# xeclass = CubicSpline(retnf[2][end:-1:begin], class_x[end:-1:begin])
# tmatclass = CubicSpline(retnf[3][end:-1:begin], class_x[end:-1:begin])
# csb2class = CubicSpline(retnf[4][end:-1:begin], class_x[end:-1:begin])
# ih = customion(𝕡, bg, xeclass, tmatclass, csb2class)

    
x_grid = bg.x_grid

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
k = 𝕡.h*kclass  #get k in our units
# k = 𝕡.h * kclass
class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))

xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))] #horizon crossing ish
println("k = ", kclass," log10k = ", log10(kclass), " h/Mpc")

#pert setup
ℓᵧ=50
ℓ_ν=50
ℓ_mν=20
reltol=1e-9
abstol=1e-9
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid)) #arrays for the massive neutrino integrated perts
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
#solve (with rsa)
results_with_rsa = boltsolve_rsa(hierarchy; reltol=reltol, abstol=abstol)
# results_with_rsa = boltsolve(hierarchy; reltol=reltol)


# Residual plot against CLASS
#reverse the arrays because ow complains
ϵ=1e-4
x_grid[1:end-1].-class_pxsnf[1,:][end:-1:1][1]
min_class_idx = minimum(findall(<(ϵ), x_grid[1:end-1]./class_pxsnf[1,:][end:-1:1][1] .- (1-ϵ) ))
x_grid_aligned = x_grid[1:end-1][min_class_idx:end]
# ours_2_Mpcm1 = 1/(2.1331e-35 *3e5) #unit conversion 1/([km/s/Mpc]*[c/km/s])


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


##
class_x = class_pxsnf[1,:][end:-1:1]

itphibolt = CubicSpline((results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]), x_grid)
itpphiclass = CubicSpline(class_pxsnf[7,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

itdeltbbolt = CubicSpline((results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]), x_grid)
itdeltbclass = CubicSpline(class_pxsnf[3,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

itpgambolt = CubicSpline(-(results_with_rsa[1,:]*4)[1:end], x_grid)
itpgamclass = CubicSpline(class_pxsnf[2,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

plot(class_x, -itpgambolt.(class_x))
plot!(class_x, -itpgamclass.(class_x), label="class")
# plot!(xlim=(-4.5,0.0), ylim=(0.9, 1.1))

##

plot()
# plot(class_x, -itdeltclass.(class_x))

class_eta = bg.η.(class_x)

plot(class_x, itphibolt.(class_x) ./ itpphiclass.(class_x), label=raw"$\Phi$")
plot!(class_x, -itdeltbbolt.(class_x) ./ itdeltbclass.(class_x), label=raw"$\delta_b$")
plot!(class_x, itpgambolt.(class_x) ./ itpgamclass.(class_x), label=raw"$\Theta_0$")

plot!(ylim=(1-1e-3, 1+1e-3), legend=:topleft, xlabel="x", ylabel="bolt / class", title="k=$(k) Mpc" * raw"$^-1$")



##

#phi
class_x = class_pxsnf[1,:][end:-1:1]
itpphiclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[7,:][end:-1:1])
plot(class_x, itpphiclass.(class_x), label="phi" )
hline!([1],ls=:dot,color=:black,label=false )


##
#matter density
# itpdelclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[4,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:])[1:end-1][min_class_idx:end]./itpdelclass.(x_grid_aligned), label="mat" )

#baryon density
# itpbarclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[3,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:])[1:end-1][min_class_idx:end]./itpbarclass.(x_grid_aligned), label="bar")
# # plot(x_grid_aligned, itpbarclass.(x_grid_aligned)./itpbarclass_rf.(x_grid_aligned), label="bar class hyrec/rf" )

#photon monopole
# itpgamclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[2,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end]./itpgamclass.(x_grid_aligned), label="pho" )

#massless neutrino monopole
# itpnu0class = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[5,:][end:-1:1])
# plot!(x_grid_aligned, -(results_with_rsa[2(ℓᵧ+1)+1,:]*4)[1:end-1][min_class_idx:end]./itpnu0class.(x_grid_aligned), label="nu0" )
# vline!([this_rsa_switch],ls=:dot,color=:green,label="rsa")
# vline!([log(1/1101)],ls=:dot,color=:green,label="~decoupling")
# vline!([log(1/3400)],ls=:dot,color=:gray,label="~MD")


vline!([xhor], color="red")
# ylims!(1.0 - 0.4e-2, 1.0 + 0.4e-2)
# xlims!(-8,0)
xlabel!(raw"$x$")
ylabel!(raw"$\frac{\delta_{i,Bolt}}{\delta_{i,CLASS}}(x)$")
title!("Compare CLASS - Bolt - k=$(@sprintf("%.3f", kclass))")
# savefig("../compare/614nonu_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

##
plot(x_grid_aligned, (results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end])
# plot!(x_grid_aligned, -itpgamclass.(x_grid_aligned), label="CLASS", ls=:dash)
plot!(class_pxsnf[1,:][end:-1:1],-class_pxsnf[2,:][end:-1:1], marker=:o)

plot!(ylim=(0.85,1.15), xlim=(-5,-4))

##

ass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[2,:][end:-1:1])
plot(bg.η.(x_grid_aligned) .* (bg.H₀*299792.458/100), -(results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end]./itpgamclass.(x_grid_aligned), 
label="pho", xscale=:log10 )

##


# plot(x_grid_aligned, itpphiclass.(x_grid_aligned), label="pho" )
# plot(x_grid_aligned, -(results_with_rsa[1,:]*4)[1:end-1][min_class_idx:end] .- itpgamclass.(x_grid_aligned), label="pho" )

##
retnf = open( @sprintf("./test/data/zack_N_class_px_nofluid_nonu_x_e.dat"),"r" ) do datafile
# an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
     [parse.(Float64, split(line)) for line in eachline(datafile)]
end

##
class_a = retnf[1]
class_x = log.(class_a)
# xeclass = LinearInterpolation(class_x[end:-1:begin],retnf[2][end:-1:begin], extrapolation_bc = Line())
# tmatclass = LinearInterpolation(class_x[end:-1:begin],retnf[3][end:-1:begin], extrapolation_bc = Line())
# csb2class = LinearInterpolation(class_x[end:-1:begin],retnf[4][end:-1:begin], extrapolation_bc = Line())


plot(Bolt.x2z.(class_x), xeclass.(class_x) ./ ih.Xₑ.(class_x), label=raw"class / bolt $X_e$")
plot!(xlim=(0,2000))
# vline!([log(1/1100)], label="recombination", xlabel="x", xlim=(-11, 0))
# plot!(class_x, )

##

plot(Bolt.x2z.(class_x), xeclass.(class_x) , label=raw"class $X_e$")
plot!(Bolt.x2z.(class_x), ih.Xₑ.(class_x), label=raw"bolt $X_e$")
vline!([log(1/1100)], label="recombination", xlabel="z", yscale=:log10, xlim=(0,2000))

##

