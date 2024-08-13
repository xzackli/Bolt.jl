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

class_x = class_pxsnf[1,:][end:-1:1]

itphibolt = CubicSpline((results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]), x_grid)
itpphiclass = CubicSpline(class_pxsnf[7,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

itdeltbbolt = CubicSpline((results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]), x_grid)
itdeltbclass = CubicSpline(class_pxsnf[3,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

itpgambolt = CubicSpline(-(results_with_rsa[1,:]*4)[1:end], x_grid)
itpgamclass = CubicSpline(class_pxsnf[2,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

plot()
class_eta = bg.η.(class_x)

plot(class_x, itphibolt.(class_x) ./ itpphiclass.(class_x), label=raw"$\Phi$")
plot!(class_x, -itdeltbbolt.(class_x) ./ itdeltbclass.(class_x), label=raw"$\delta_b$")
plot!(class_x, itpgambolt.(class_x) ./ itpgamclass.(class_x), label=raw"$\Theta_0$")

plot!(ylim=(1-1e-3, 1+1e-3), legend=:topleft, xlabel="x", ylabel="bolt / class", title="k=$(k) Mpc" * raw"$^-1$")

