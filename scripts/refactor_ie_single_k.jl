using Bolt
using Plots
using Plots.PlotMeasures
using NumericalIntegration
using Printf
using DelimitedFiles
using Interpolations
using BenchmarkTools

using Bolt: spline #FIXME why do I have to import this here but NOT in bg?

# Load some saved hierarchy answers to compare against (and start from)
k_options = ["p03", "p3", "1p0"] #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
k_choice = k_options[2]
kMpc = parse(Float64, replace(k_choice,"p"=>".")) #DUMB does not work for comparison
ret = readdlm( @sprintf("./test/data/bolt_px_k%s_fine.dat",k_choice) )
retnf_class = open( @sprintf("./test/data/class_px_k%s_nofluid_re.dat",k_choice),"r" ) do datafile
# an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#the second column is just a repeated k value, so remember it and delete col
kclass = retnf_class[2][1] #read class k mode from file (in h/Mpc)
dx = ret[2,1]-ret[1,1]

#generate some background/ionization history
𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q);
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b);
ih = IonizationHistory(𝕣, 𝕡, bg);
Mpcfac = bg.H₀*299792.458/100.
k = Mpcfac*kclass #get k in our units

#input to the ie integrator struct (akin to hierarchy)
ℓᵧ=2
ℓ_ν=50
ℓ_mν=20
reltol=1e-6 #cheaper  rtol

# Initial setup
N_iters = 5 # in principle should replace this with a tolerance criterion
x₀ = bg.x_grid
Θ₂₀,Π₀ = zeros(length(bg.x_grid)),zeros(length(bg.x_grid)) #SHOULD use something better
ie_0 = IE(BasicNewtonian(), 𝕡, bg, ih, k,
        linear_interpolation(x₀,Θ₂₀),linear_interpolation(x₀,Θ₂₀),
        20,100,200,
        ℓ_ν, ℓ_mν, n_q);
u_ie = itersolve(N_iters,ie_0);

plot(ret[:,1],ret[:,1+1],label="hierarchy")
plot!(x_grid_ie(ie_0),u_ie[1,:],ls=:dash,label="iter 5")


