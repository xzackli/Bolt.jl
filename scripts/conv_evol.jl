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


𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k = (bg.H₀*299792.458/100)*kclass #get k in our units
ℓᵧ=50
ℓ_ν=50
ℓ_mν=20
reltol=1e-8 #cheaper  rtol

# x-evolution
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, 50, ℓ_ν, ℓ_mν,n_q)
results = boltsolve(hierarchy;reltol=reltol)

# η-evolution
η2x = linear_interpolation(bg.η,bg.x_grid)
hierarchy_conf = ConformalHierarchy(hierarchy,η2x)
results_conf = boltsolve_conformal(hierarchy_conf;reltol=reltol)

plot(bg.η.(results.t)*bg.H₀*3e5/100,results[1,:],label="hierarchy-x-1e-8",color=:blue)
plot!(results_conf.t,results_conf[1,:],label="hierarchy-η-1e-8",color=:orange,ls=:dash,xscale=:log10)
plot!(bg.η.(retnf_class[1,:])*bg.H₀*3e5/100,-retnf_class[2+1,:]/4/𝕡.h,label="hierarchy-CLASS",color=:green,ls=:dot)

xlims!(bg.η(-12)*bg.H₀*3e5/100,bg.η(0)*bg.H₀*3e5/100)
xlabel!("η")
ylabel!("temp mono")
savefig("../temp_mono_hier_xeta_rtol-8.png")
