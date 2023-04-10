# Dev script for combining photons and neutrinos
using Bolt
using Plots
using Plots.PlotMeasures
using NumericalIntegration
using Printf
using DelimitedFiles
using Interpolations
using BenchmarkTools


# /// Setup ///
k_options = ["p03", "p3", "1p0"] #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
k_choice = k_options[1]
kMpc = parse(Float64, replace(k_choice,"p"=>".")); 
ret = readdlm( @sprintf("./test/data/bolt_px_k%s_fine.dat",k_choice) );
retnf_class = open( @sprintf("./test/data/class_px_k%s_nofluid_re.dat",k_choice),"r" ) do datafile
# an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#the second column is just a repeated k value, so remember it and delete col
kclass = retnf_class[2][1] #read class k mode from file (in h/Mpc)
dx = ret[2,1]-ret[1,1]
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))


𝕡 = CosmoParams();
# 𝕡 = CosmoParams(
#     h = 0.6774,  # hubble factor
#     Ω_b = 0.0486, 
#     Ω_m = 0.2589,
#     Σm_ν = 0.15
# ) # Planck15 modifications to h, Ω_b,Ω_c, make mnu=0 
n_q=15
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q);
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b); #FIXME γΩ
ih = IonizationHistory(𝕣, 𝕡, bg);
Mpcfac = bg.H₀*299792.458/100.
k = Mpcfac*kclass #get k in our units
# Hierarchy for comparison purposes
ℓᵧ=50
ℓ_mν=20
ℓ_ν=50
pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
reltol=1e-12
#solve the hierarchy just to be sure
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q);
results=zeros(pertlen,length(bg.x_grid));
perturb = boltsolve(hierarchy;reltol=reltol);

for (i_x, x) in enumerate(bg.x_grid)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
end
#conformal hierarchy
η2x = linear_interpolation(bg.η,bg.x_grid);
hierarchy_conf = ConformalHierarchy(hierarchy,η2x);
perturb_conf = boltsolve_conformal(hierarchy_conf;reltol=reltol);
results_conf=zeros(pertlen,length(perturb_conf.t));
for (i_t, t) in enumerate(perturb_conf.t)
    u = perturb_conf(t)  #z this can be optimized away, save timesteps at the grid!
    results_conf[:,i_t] = u #z should use unpack somehow
end

# writedlm("./test/data/Bolt_allperts_nonu_pholmax$(ℓᵧ)_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat",
#           hcat(η2x(perturb_conf.t),
#           results_conf[3,:], results_conf[3,:] .+ results_conf[(ℓᵧ+1)+1,:] .+ results_conf[(ℓᵧ+1)+3,:],
#           results_conf[2(ℓᵧ+1)+1,:],results_conf[2(ℓᵧ+1)+3,:],
#           [results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:] for idx_q in 1:n_q]...,
#           [results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+idx_q,:] for idx_q in 1:n_q]...))



planck_heavynu_ansatz_data = readdlm("./test/data/Bolt_allperts_mnu0p15_pholmax$(ℓᵧ)_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat");



planck_pho_heavynu_ansatz_Θ₂ = linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,2]);
planck_pho_heavynu_ansatz_Π = linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,3]);
planck_nu_heavynu_ansatz₀ = [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,4]),
                        [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,1+4+idx_q]) for idx_q in 1:n_q]...];

planck_nu_heavynu_ansatz₂ = [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,5]),
                        [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,1+4+n_q+idx_q]) for idx_q in 1:n_q]...];


# Initial setup
N_iters = 2 # in principle should replace this with a tolerance criterion

# Θ₂_interp,Π₂_interp =  linear_interpolation(η2x(perturb_conf.t),results_conf[3,:]),linear_interpolation(η2x(perturb_conf.t),results_conf[3,:]+results_conf[(ℓᵧ+1)+1,:]+results_conf[(ℓᵧ+1)+3,:])
# all_interps₀ = [linear_interpolation(η2x(perturb_conf.t),results_conf[2(ℓᵧ+1)+1,:]),[linear_interpolation(η2x(perturb_conf.t),results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]) for idx_q in 1:n_q]...];
# all_interps₂ = [linear_interpolation(η2x(perturb_conf.t),results_conf[2(ℓᵧ+1)+3,:]),[linear_interpolation(η2x(perturb_conf.t),results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+idx_q,:]) for idx_q in 1:n_q]...];


ie_0 = IEγν(BasicNewtonian(), 𝕡, bg, ih, k,
        # Θ₂_interp,Π₂_interp,
        # all_interps₀,all_interps₂,
        planck_pho_heavynu_ansatz_Θ₂,planck_pho_heavynu_ansatz_Π,
        planck_nu_heavynu_ansatz₀,planck_nu_heavynu_ansatz₂,
        300,300,800, #exact optimal choice of these is k,η-dependent...
        n_q);

η_switch =1.0;
switch_idx = argmin(abs.(bg.η .-η_switch))
η2x_late = linear_interpolation(bg.η.(bg.x_grid[switch_idx:end]), bg.x_grid[switch_idx:end]);
x2η_late = linear_interpolation( bg.x_grid[switch_idx:end],bg.η.(bg.x_grid[switch_idx:end]));
cie_0 = ConformalIEγν(ie_0,η2x_late);
η_switch_use = bg.η[switch_idx];
u0_ie = get_switch_u0(η_switch_use,hierarchy_conf,reltol);
M = 2048*4


h_boltsolve_conformal_flex(hierarchy_conf, η_switch_use, hierarchy_conf.hierarchy.bg.η(hierarchy_conf.hierarchy.bg.x_grid[end]),  
                            initial_conditions(hierarchy_conf.hierarchy.bg.x_grid[1], hierarchy_conf.hierarchy),reltol=reltol);


xx_k,Θ₂,Π,𝒳₀_k,𝒳₂_k,perturb_k = itersolve(N_iters,cie_0,M,bg.x_grid[switch_idx],0.0,u0_ie;reltol=reltol);

#Accuracy/timing

@btime get_switch_u0(η_switch_use,hierarchy_conf,reltol);
#1.017 s (19372055 allocations: 439.42 MiB)
#Planck Ansatz (PA), Niter=2:    675.599 ms (19372055 allocations: 439.42 MiB)
@btime itersolve(N_iters,cie_0,M,bg.x_grid[switch_idx],0.0,u0_ie;reltol=reltol);
#11.541 s (78299871 allocations: 2.57 GiB)
#PA,Niter=2:   3.369 s (36001249 allocations: 1.14 GiB) #everything should be lower here by a factor of 2/5, which is about right...
@btime boltsolve_conformal(hierarchy_conf;reltol=reltol);
#4.814 s (86049198 allocations: 1.73 GiB)
#PA,Niter=2 (same code as above but running at simlar time):  3.394 s (86049198 allocations: 1.73 GiB)


#flamegraphs
@profview get_switch_u0(η_switch_use,hierarchy_conf,reltol);
#most time in in the Jacobian as expected for hierarchy solve, lowest level time is dominated by hierarchy derivative computations themselves
@profview itersolve(N_iters,cie_0,M,bg.x_grid[switch_idx],0.0,u0_ie;reltol=reltol);
#about half the time in boltsolve_conformal_flex, 1/3 in photon part, 1/6 in neutrino part
# within these
# 1) In nlsolve (which I guess is for jac?), 1/3 of the ie_conformal time is in ρ_σ, rest is just elsewhere in ie!
# 2) All in g_weight_trapz_ie, 2/3 in _IΘ2 (mostly R2 and bspline), 1/4 in _IΠ, rest in W±
# 2) All in g_weight_trapz_ie, 2/3 in _IΘ2 (mostly R2 and bspline), 1/4 in _IΠ, rest in W±
# 3) Some easy optimizations with "materialize", ~1/4 in χνz, rest of most in FFTs
@profview boltsolve_conformal(hierarchy_conf;reltol=reltol);
# again nlsolve jacobian dominates, and at bottom have hierarchy cost, looks like there is some broadcasting here that might be eliminated through

# use monopole l2 error as an accuracy diagnostic
sqrt(sum( (perturb_k(bg.η(xx_k))[1,:].- perturb_conf(bg.η(xx_k))[1,:]).^2 ))/length(xx_k)

#Checks:

plot(perturb_conf.t,results_conf[3,:],color=:black,xscale=:log10)
plot!(bg.η(xx_k),Θ₂(xx_k),ls=:dash)
xlims!(η_switch,bg.η[end])
ylabel!("Θ₂(η)")
xlabel!("η")






plot(bg.η(bg.x_grid[2:end]),ih.τ(bg.x_grid[1:end-1]) .-ih.τ(bg.x_grid[2:end]) .+1e-16,xscale=:log10,yscale=:log10)
hline!([1e-3])
#start with massless case
#back eta
bη =  bg.η[end] .- bg.η  
plot!(bg.η(bg.x_grid[2:end]),(bη[1:end-1] .-bη[2:end] ) .+1e-16,xscale=:log10,yscale=:log10)
# plot!(bg.η(bg.x_grid[2:end]),-(1.0./bη[1:end-1] .- 1.0./bη[2:end] ) .+1e-16,xscale=:log10,yscale=:log10)
plot!(bg.η(bg.x_grid[2:end]),(1.0./(bη[1:end-1] .-bη[2:end]) ) .+1e-16,xscale=:log10,yscale=:log10)

Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
logqmin,logqmax=log10(Tν/30),log10(Tν*30)
q15 = Bolt.xq2q(bg.quad_pts[15],logqmin,logqmax)
χνs15 = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q15,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts,𝕡) for x in bg.x_grid])*Mpcfac
q1 = Bolt.xq2q(bg.quad_pts[1],logqmin,logqmax)
χνs1 = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q1,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts,𝕡) for x in bg.x_grid])*Mpcfac
plot(bg.η,χνs .+1e-16,xscale=:log10,yscale=:log10)
bχνs15 = χνs15[end] .- χνs15
bχνs1 = χνs1[end] .- χνs1

plot!(bg.η(bg.x_grid[2:end]),1.0./(bχνs15[1:end-1] .-bχνs15[2:end] ) .+1e-16,xscale=:log10,yscale=:log10)
plot!(bg.η(bg.x_grid[2:end]),1.0./(bχνs1[1:end-1] .-bχνs1[2:end] ) .+1e-16,xscale=:log10,yscale=:log10)

plot(bg.η, Bolt.j0.(bχνs15/Mpcfac*k))
plot(bg.η, abs.(Bolt.j0.(bχνs15/Mpcfac*k)),yscale=:log10)
plot!(bg.η(bg.x_grid[2:end]),ih.τ(bg.x_grid[1:end-1]) .-ih.τ(bg.x_grid[2:end]) .+1e-16,xscale=:log10,yscale=:log10)

Bolt.j0.(bχνs15/Mpcfac*k)
bg.η(bg.x_grid[2:end])

SinInteg