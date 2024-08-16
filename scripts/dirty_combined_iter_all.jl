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
k_choice = k_options[3]
kMpc = parse(Float64, replace(k_choice,"p"=>".")); #DUMB does not work for comparison
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
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3)/2.0:ret[end,1], nq=n_q);
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
perturb = boltsolve(hierarchy; abstol=1e-6,reltol=reltol);
for (i_x, x) in enumerate(bg.x_grid)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
end
#conformal hierarchy
η2x = linear_interpolation(bg.η,bg.x_grid);
hierarchy_conf = ConformalHierarchy(hierarchy,η2x);
perturb_conf = boltsolve_conformal(hierarchy_conf;reltol=reltol,abstol=1e-6);
results_conf=zeros(pertlen,length(perturb_conf.t));
for (i_t, t) in enumerate(perturb_conf.t)
    u = perturb_conf(t)  #z this can be optimized away, save timesteps at the grid!
    results_conf[:,i_t] = u #z should use unpack somehow
end

plot(bg.η , results[3,:],ls=:dash)
plot!(perturb_conf.t,results_conf[3,:])
xlims!(0.0,500.0)

plot(bg.η , results[1,:],ls=:dash,xscale=:log10)
plot!(perturb_conf.t,results_conf[1,:])
xlims!(.1,100.0)


plot(bg.η , results[end-3,:],ls=:dash,xscale=:log10)
plot!(perturb_conf.t,results_conf[1,:])
xlims!(.1,100.0)


@time boltsolve(hierarchy; reltol=reltol);
@btime boltsolve(hierarchy; reltol=reltol);
#@btime 3.902 s (114143281 allocations: 2.47 GiB)
@btime boltsolve_conformal(hierarchy_conf;reltol=reltol);
#@btime 4.931 s (86047329 allocations: 1.72 GiB)


function get_perts(u,ie::IEγν{T},x) where T
    k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie) 
    Θ₂ = ie.sΘ2(x)
    𝒩₀,𝒩₂  = ie.s𝒳₀[1](x),ie.s𝒳₂[1](x)
    ℳ₀,ℳ₂ = zeros(T,nq),zeros(T,nq)
    for idx_q in 1:nq
        ℳ₀[idx_q] = ie.s𝒳₀[idx_q+1](x)
        ℳ₂[idx_q] = ie.s𝒳₂[idx_q+1](x)
    end
    ρℳ, σℳ  =  @views ρ_σ(ℳ₀, ℳ₂, bg, a, par)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ₂ +
                                  Ω_ν * 𝒩₂
                                  + σℳ / bg.ρ_crit /4
                                  )
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩₀
        + a^(-2) * ρℳ / bg.ρ_crit
        )

    Π = ie.sΠ(x)
    return Φ′,Ψ,Θ[0],Π,v_b
end

#merged iterate
function iterate(Θ₂_km1,Π_km1, 𝒳₀_km1,𝒳₂_km1, 
        ie::IEγν{T},
        #  𝕡::CosmoParams{T}, bg, ih, k, n_q,
        M,x_ini, x_fin,u0,reltol) where T
    𝕡, bg, ih, k, n_q,Nᵧ₁,Nᵧ₂,Nᵧ₃ = ie.par,ie.bg,ie.ih,ie.k,ie.nq,ie.Nᵧ₁,ie.Nᵧ₂,ie.Nᵧ₃ #FIXME get rid of this line
    Mpcfac = bg.H₀*299792.458/100.
    ie_k = IEγν(BasicNewtonian(), 𝕡, bg, ih, k,
            Θ₂_km1,Π_km1,
            𝒳₀_km1,𝒳₂_km1,
            Nᵧ₁,Nᵧ₂,Nᵧ₃,
            n_q)

    # Do the truncated boltzmann solve
    perturb_k = boltsolve_flex(ie_k, x_ini, x_fin, u0; reltol=reltol)

    # Get metric and photon-relevant perturbation variables
    #FIXME may want to put this in its own function
    xgi = x_grid_ie(ie_k,x_ini,x_fin)
    
    N = length(xgi)
    Φ′,Ψ,Θ₀,Π,v_b = zeros(N),zeros(N),zeros(N),zeros(N),zeros(N)
    # for (j,u) in enumerate( eachcol(u_all_k) )
    u_all = perturb_k(xgi)

    for (j,u) in enumerate( eachcol(u_all) )
        Φ′[j],Ψ[j],Θ₀[j],Π[j],v_b[j] = get_perts(u,ie_k,xgi[j])
    end

    #photons
    aΘ₂_k,aΠ_k = zeros(N),zeros(N)
    for i in 3:N
        aΘ₂_k[i],aΠ_k[i] = g_weight_trapz_ie(i,xgi,ie_k,Φ′,Ψ,Θ₀,Π,v_b)
    end
    Θ₂_k,Π_k = linear_interpolation(xgi,aΘ₂_k), linear_interpolation(xgi,aΠ_k)
    sΦ′,sΨ = linear_interpolation(xgi,Φ′),linear_interpolation(xgi,Ψ)

    #neutrinos
    xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k,M,u0,perturb_k.t,sΦ′,sΨ) 
    return xx,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k
end
function iterate(Θ₂_km1,Π_km1, 𝒳₀_km1,𝒳₂_km1, 
                #  ie::IEγν{T},
                cie::ConformalIEγν{T},
                #  𝕡::CosmoParams{T}, bg, ih, k, n_q,
                 M,x_ini, x_fin,u0,reltol) where T
    ie = cie.ie
    𝕡, bg, ih, k, n_q,Nᵧ₁,Nᵧ₂,Nᵧ₃ = ie.par,ie.bg,ie.ih,ie.k,ie.nq,ie.Nᵧ₁,ie.Nᵧ₂,ie.Nᵧ₃ #FIXME get rid of this line
    Mpcfac = bg.H₀*299792.458/100.
    # println("photon types: ", typeof(Θ₂_km1),typeof(Π_km1))
    # println("neutrino types: ", typeof(𝒳₀_km1),typeof(𝒳₂_km1))
    # println("2: ", 𝒳₀_km1[1](-20.0))
    ie_k = IEγν(BasicNewtonian(), 𝕡, bg, ih, k,
                     Θ₂_km1,Π_km1,
                     𝒳₀_km1,𝒳₂_km1,
                     Nᵧ₁,Nᵧ₂,Nᵧ₃,
                     n_q)

    cie_k = ConformalIEγν(ie_k,cie.η2x)
    # println("3: ", 𝒳₀_km1[1](-20.0))
    # println("4: ", ie_k.s𝒳₀[1](-20.0))

    # println("xini: ", ie_k.sΘ2(x_ini))

    # Do the truncated boltzmann solve
    # perturb_k = boltsolve_flex(ie_k, x_ini, x_fin, u0; reltol=reltol)
    # println("ie.sΘ2(xini) = ",ie_k.sΘ2(x_ini))
    # println("ie.sΘ2(xini) = ",ie_k.sΘ2(cie_k.η2x(bg.η(x_ini))))


    #need to fix this first...
    # println("flop death operation: ", bg.η(x_fin), ", ", (bg.η(x_fin)*Mpcfac)  / Mpcfac )
    # println("pre-death operation: ", cie_k.η2x( bg.η(x_fin)  ) )
    # println("death operation: ", cie_k.η2x( ( bg.η(x_fin)*Mpcfac)  / Mpcfac ) )

    perturb_k = boltsolve_conformal_flex(cie_k, bg.η(x_ini), bg.η(x_fin), u0, reltol=reltol)

    # cie_k = ConformalIEγν(ie_k,linear_interpolation(perturb_k.t/Mpcfac,cie.η2x(perturb_k.t/Mpcfac)))
    # x_ini = cie.η2x(perturb_k.t[1]/Mpcfac)

    # println("cie.η2x(perturb_k.t[1] ./Mpcfac)", cie.η2x(perturb_k.t[1] /Mpcfac))
    # println("cie.η2x.(perturb_k.t ./Mpcfac)[1]", cie.η2x.(perturb_k.t ./Mpcfac)[1])
    # println("cie.η2x(perturb_k.t ./Mpcfac)[1]",cie.η2x(perturb_k.t ./Mpcfac)[1])
    # println("cie.η2x(bg.η(x_ini))",cie.η2x(bg.η(x_ini)))
    # println("etas: ", bg.η(x_ini)*Mpcfac, ", ",perturb_k.t[1])
    #this is drivig me crazy - iterpolatio is not closed...
    # FIXME The η2x(perturb_k.t[1]) is not the same as x_ini for any choice of the above
    
    # Get metric and photon-relevant perturbation variables
    #FIXME may want to put this in its own function
    
    
    # xgi = x_grid_ie(ie_k,x_ini,x_fin)
    xgi = cie_k.η2x(η_grid_ie(ie_k,bg.η(x_ini),bg.η(x_fin),cie_k.η2x)) #FIXME HACK


    # println("xgi[1], xgi[end]",xgi[1],xgi[end])
    # println("x_ini,x_fin",x_ini,", ",x_fin)
    N = length(xgi)
    # assert(length(xgi)==length(Π_k))
    Φ′,Ψ,Θ₀,Π,v_b = zeros(N),zeros(N),zeros(N),zeros(N),zeros(N)
    # for (j,u) in enumerate( eachcol(u_all_k) )
    #why is perturb_k allowed to be longer than N? probably because I don't enforce save_at_grid...
    # println("ptk size: ",size(perturb_k))
    # u_all = perturb_k(xgi)
    # u_all = perturb_k(bg.η(xgi)*Mpcfac)
    u_all = perturb_k(bg.η(xgi))

    # println("u_all size: ",size(u_all))
    for (j,u) in enumerate( eachcol(u_all) )
            Φ′[j],Ψ[j],Θ₀[j],Π[j],v_b[j] = get_perts(u,ie_k,xgi[j])
    end

    #photons
    aΘ₂_k,aΠ_k = zeros(N),zeros(N)
    for i in 3:N
            aΘ₂_k[i],aΠ_k[i] = g_weight_trapz_ie(i,xgi,ie_k,Φ′,Ψ,Θ₀,Π,v_b)
    end
    Θ₂_k,Π_k = linear_interpolation(xgi,aΘ₂_k), linear_interpolation(xgi,aΠ_k)
    sΦ′,sΨ = linear_interpolation(xgi,Φ′),linear_interpolation(xgi,Ψ)
    
    #neutrinos
    # xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k,M,u0,perturb_k.t,sΦ′,sΨ) 
    # println("xgi[1]", xgi[1])
    # println("cie.η2x(perturb_k.t ./Mpcfac)[1]",cie.η2x(perturb_k.t ./Mpcfac)[1])
    # xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k,M,u0,cie.η2x(perturb_k.t ./Mpcfac),sΦ′,sΨ) 
    xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k,M,u0,xgi,sΦ′,sΨ)

    return xx,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k
end

# Merged itersolve
function itersolve(Nₖ::Int,
                    ie_0::IEγν{T},
                    M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
                ie_0 = cie_0.ie
                𝒳₀_k,𝒳₂_k = ie_0.s𝒳₀,ie_0.s𝒳₂
                Θ₂_k,Π_k = ie_0.sΘ2,ie_0.sΠ
                xx_k,perturb_k = nothing,nothing
                for k in 1:Nₖ
                println("iter = ",k)
                xx_k,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate(Θ₂_k,Π_k, 𝒳₀_k,𝒳₂_k,
                                                        ie_0,
                                                        M,x_ini,x_fin,u0,
                                                        reltol)
                end
                return xx_k,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k
end

function itersolve(Nₖ::Int,
                #    ie_0::IEγν{T},
                   cie_0::ConformalIEγν{T},
                   M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
    ie_0 = cie_0.ie
    𝒳₀_k,𝒳₂_k = ie_0.s𝒳₀,ie_0.s𝒳₂
    # println("1: ", 𝒳₀_k[1](-20.0))
    Θ₂_k,Π_k = ie_0.sΘ2,ie_0.sΠ
    xx_k,perturb_k = nothing,nothing
    for k in 1:Nₖ
        println("iter = ",k)
        xx_k,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate(Θ₂_k,Π_k, 𝒳₀_k,𝒳₂_k,
                                                    cie_0,# ie_0,
                                                    M,x_ini,x_fin,u0,
                                                    reltol)
    end
    return xx_k,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k
end


# Helper functon for switch
# function get_switch_u0(η,hierarchy,reltol) #Input is η of the switch
function get_switch_u0(η,hierarchy_conf,reltol) 
    # This function assumes truncated hierarchies for all neutrinos (but not yet photons)
    hierarchy = hierarchy_conf.hierarchy
    bg =hierarchy.bg
    Mpcfac = bg.H₀*299792.458/100.
    # switch_idx = argmin(abs.(bg.η*Mpcfac .-η)) #for now we use the bg to find the switch
    switch_idx = argmin(abs.(bg.η .-η)) #for now we use the bg to find the switch
    #solve the split ode
    ℓᵧ,ℓ_ν,n_q = hierarchy.ℓᵧ,hierarchy.ℓ_ν, hierarchy.nq
    pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
    # \/ we want to report this timing to get a full picture of total time (early+late)
    # sol_early_c = h_boltsolve_flex(hierarchy, bg.x_grid[1], bg.x_grid[switch_idx],  initial_conditions(bg.x_grid[1], hierarchy),reltol=reltol);
    sol_early_c = Bolt.h_boltsolve_conformal_flex(hierarchy_conf, bg.η(bg.x_grid[1]), bg.η(bg.x_grid[switch_idx]),  initial_conditions(bg.x_grid[1], hierarchy),reltol=reltol);

    # Get the new initial conditions
    u0_ie = zeros(2(2) + (0+1) + (0+1)*n_q + 5);
    # The first split for photons
    u0_ie[1] = sol_early_c.u[end][1]
    u0_ie[2] = sol_early_c.u[end][2]
    u0_ie[3] = sol_early_c.u[end][(ℓᵧ+1)+1]
    u0_ie[4] = sol_early_c.u[end][(ℓᵧ+1)+3]
    #set the massless neutrino dipole
    u0_ie[2(2)+1] = sol_early_c.u[end][2(ℓᵧ+1)+2]

    #massive neutrinos, now we just do the dipole again
    # start at the dipole first q idx, go up to the last dipole q idx (in the hierarchy)   
    for i in 1:n_q 
        u0_ie[2(2)+1+i] = sol_early_c.u[end][2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1+i]
    end
    #metric and cold perts
    for i in 1:5 #skip the higher massless hierarchy multipoles
        u0_ie[2(2)+1+n_q+i] = sol_early_c.u[end][pertlen-5+i]
    end
    return u0_ie
end

# η2x(2π/k/300)
# 2π/k/300*Mpcfac
# No cosmological modes of interest (up to k of 10) 
#enter the horizon before x=-13.3? which is 2π/k/300, ctime ~0.7 Mpc?
#^So for any reasonable k we should start before "horizon entry" so defined


# Initial setup
N_iters = 3 # in principle should replace this with a tolerance criterion
#sanity check (full hierarchy ansatz)
x_Θ₂_interp,x_Π₂_interp =  linear_interpolation(collect(bg.x_grid),results[3,:]),linear_interpolation(collect(bg.x_grid),results[3,:]+results[(ℓᵧ+1)+1,:]+results[(ℓᵧ+1)+3,:])
x_all_interps₀ = [linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+1,:]),[linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]) for idx_q in 1:n_q]...];
x_all_interps₂ = [linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+3,:]),[linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+idx_q,:]) for idx_q in 1:n_q]...];

ie_0 = IEγν(BasicNewtonian(), 𝕡, bg, ih, k,
        x_Θ₂_interp,x_Π₂_interp,
        x_all_interps₀,x_all_interps₂,
        300,300,800,
        n_q);
# constructor at least works initially, good...
η_switch =1.0;
switch_idx = argmin(abs.(bg.η .-η_switch))
η2x_late = linear_interpolation(bg.η.(bg.x_grid[switch_idx:end]), bg.x_grid[switch_idx:end])
x2η_late = linear_interpolation( bg.x_grid[switch_idx:end],bg.η.(bg.x_grid[switch_idx:end]))
cie_0 = ConformalIEγν(ie_0,η2x_late);


bg.x_grid[switch_idx:end]
bg.η.(bg.x_grid[switch_idx:end])
x2η_late(η2x_late(bg.η(bg.x_grid[switch_idx]))) == bg.η(bg.x_grid[switch_idx])
bg.η(η2x_late(bg.η.(bg.x_grid[switch_idx]))) == bg.η(bg.x_grid[switch_idx])


#make it so we actually hit zero...
# xgi_test = x_grid_ie(ie_0,η2x(1.0/Mpcfac),0.0)
xgi_test = x_grid_ie(ie_0,η2x(1.0),0.0)
ηgi_test = η_grid_ie(ie_0,η_switch_use,bg.η[end],η2x)
ηgi_test
xgi_test
η2x(ηgi_test)

bg.η[end]
2.2619502561780378e33
η2x_late(    2.2619502561780378e33)
η2x_late(    bg.η[end])
η2x_late(    bg.η(bg.x_grid[end]))
bg.η(bg.x_grid[end]) == bg.η[end]


plot(bg.x_grid,abs.(ih.g̃).+1e-16,yscale=:log10)
vline!([log(1/(1300+1)),log(1/(800+1))])
ylims!(1e-5,1e1)

plot(bg.x_grid,abs.(ih.τ).+1e-16,yscale=:log10)
vline!([log(1/(1300+1)),log(1/(800+1))])



η2x_late(2.2619502561780378e33)
η2x_late(bg.η[end])
η2x_late(bg.η(bg.x_grid[end]))

η_switch_use = bg.η[switch_idx];#1.0;
u0_ie = get_switch_u0(η_switch_use,hierarchy_conf,reltol);
M = 2048*4

xx_kx,Θ₂x,Πx,𝒳₀_kx,𝒳₂_kx,perturb_kx = itersolve(N_iters,ie_0,M,bg.x_grid[switch_idx],-5.0,u0_ie;reltol=reltol);
xx_k,Θ₂,Π,𝒳₀_k,𝒳₂_k,perturb_k = itersolve(N_iters,cie_0,M,bg.x_grid[switch_idx],-5.0,u0_ie;reltol=reltol);

ηlateunit==bg.η[switch_idx]
ηlateunit ≈ bg.η[switch_idx]
ηlateunit, bg.η[switch_idx]

#how is this possible? This is now on a node of η2x, so it can't be wrong??
xinit=η2x_late()
ηlateunit = η_switch/Mpcfac
η2x_late[1]
η2x_late(bg.η(xinit))==xinit
η2x_late(bg.η(xinit)) ≈ xinit
η2x_late(bg.η(xinit)),xinit

# does this at least close? no??
x2η_late(η2x_late(ηlateunit)) == ηlateunit
x2η_late(η2x_late(ηlateunit)) ≈ ηlateunit
x2η_late(η2x_late(ηlateunit)), ηlateunit

# full reltol (1e-12)
@btime h_boltsolve_flex(hierarchy, bg.x_grid[1], bg.x_grid[argmin(abs.(bg.η*Mpcfac .-η_switch))],  initial_conditions(bg.x_grid[1], hierarchy),reltol=reltol);
#@btime  2.760 s (56284055 allocations: 1.24 GiB)
@btime h_boltsolve_conformal_flex(hierarchy_conf, bg.η[1], η_switch/Mpcfac,  initial_conditions(bg.x_grid[1], hierarchy),reltol=reltol);
#@btime  973.598 ms (19371946 allocations: 439.20 MiB)

#cheap reltol
@btime h_boltsolve_flex(hierarchy, bg.x_grid[1], bg.x_grid[argmin(abs.(bg.η*Mpcfac .-η_switch))],  initial_conditions(bg.x_grid[1], hierarchy),reltol=1e-6);
#@btime 2.688 s (55033855 allocations: 1.21 GiB)
@btime h_boltsolve_conformal_flex(hierarchy_conf, bg.η[1], η_switch/Mpcfac,  initial_conditions(bg.x_grid[1], hierarchy),reltol=1e-6);
#@btime  826.294 ms (16888485 allocations: 383.91 MiB)

@btime get_switch_u0(η_switch,hierarchy,reltol);
#@btime 2.678 s (56284088 allocations: 1.24 GiB) #1.0 switch
#@btime 3.238 s (64638594 allocations: 1.42 GiB) #5.0 switch
@btime itersolve(N_iters,ie_0,M,bg.x_grid[switch_idx],0.0,u0_ie;reltol=reltol);
#@btime 4.000 s (28954162 allocations: 1.34 GiB) #1.0 switch
#@btime 2.887 s (25420050 allocations: 1.25 GiB) #5.0 switch
# Looks like the getswitch would be a factor of maybe 3 faster if we used conformal solve...

plot!(xx_k,Θ₂(xx_k),label="IE-η")
plot!(xx_k,Θ₂x(xx_k),label="IE-x",ls=:dash)
plot(xx_k,ie_0.sΘ2(xx_k),label="hierarchy",color=:black)
vline!([log(1/(1300+1)),log(1/(800+1))])
xlabel!("ln(a)")
ylabel!("Θ₂")
xlims!(-10,-6)
       
plot!(bg.η(xx_k),Θ₂(xx_k),ls=:dash)
plot!(bg.η(xx_k),Θ₂x(xx_k),ls=:dash)
plot!(bg.η(xx_k),ie_0.sΘ2(xx_k))
xlims!(1,450)

plot(xx_k,Π(xx_k))
plot!(xx_kx,Πx(xx_kx))
plot!(xx_k,ie_0.sΠ(xx_k))

plot(xx_k,𝒳₀_k[1](xx_k))
plot!(xx_kx,𝒳₀_kx[1](xx_kx))
plot!(xx_k,ie_0.s𝒳₀[1](xx_k))
ylims!(-0.1,0.1)
plot(xx_k,𝒳₂_k[1](xx_k))
plot!(xx_kx,𝒳₂_kx[1](xx_kx))
plot!(xx_k,ie_0.s𝒳₂[1](xx_k))
plot(xx_k,abs.(𝒳₀_k[2](xx_k)),yscale=:log10)
plot!(xx_kx,abs.(𝒳₀_kx[2](xx_kx)))
plot!(xx_k,abs.(ie_0.s𝒳₀[2](xx_k)))
plot(xx_k,abs.(𝒳₂_k[2](xx_k)),yscale=:log10)
plot!(xx_kx,abs.(𝒳₂_kx[2](xx_kx)))
plot!(xx_k,abs.(ie_0.s𝒳₂[2](xx_k)))

# plot(xx_k,perturb_k(xx_k)[1,:])
plot(xx_k,perturb_k(bg.η(xx_k))[1,:],label="hr η-photon grid")#    ls=:dash)
plot!(bg.η(xx_k),perturb_k(bg.η(xx_k))[1,:],label="η-photon grid")#,xscale=:log10)#    ls=:dash)
plot!(bg.η,results[1,:])
plot(xx_k,perturb_k(bg.η(xx_k))[1,:],label="x-photon grid",ls=:dash)
plot!(xx_kx,perturb_kx(xx_kx)[1,:],label="x-photon grid",ls=:dash)
plot!(bg.η(xx_kx),perturb_kx(xx_kx)[1,:],label="x-photon grid",ls=:dash)

plot!(bg.x_grid,results[1,:])
# plot!(η2x(perturb_conf.t),results_conf[1,:],label="conf hierarchy")
plot!(perturb_conf.t,results_conf[1,:],label="conf hierarchy")
xlims!(-13.0,-6.0)

xlims!(5000,8000)

xlims!(1,500)

ylims!(-.25,.25)

vline!([log10(1/(1300+1)),log10(1/(800+1))])


xhor = bg.x_grid[argmin(abs.(k .* bg.η /Mpcfac .- 2π))]
xdec = bg.x_grid[argmin(abs.( -ih.τ′ .* bg.ℋ .*bg.η /Mpcfac .- 1))]

cie_0.η2x(η_grid_ie(ie_0,bg.η[switch_idx],bg.η(-5.0),cie_0.η2x)) == sort(cie_0.η2x(η_grid_ie(ie_0,bg.η[switch_idx],bg.η(-5.0),cie_0.η2x)))

#MB version of this plot
plot(xx_k/log(10),abs.(perturb_k(bg.η(xx_k))[1,:]),label="η-photon grid",yscale=:log10)
plot!(xx_k/log(10),abs.(perturb_kx(xx_kx)[1,:]),label="x-photon grid",ls=:dash)
plot!( η2x(perturb_conf.t)/log(10),abs.(results_conf[1,:]),label="conf hierarchy")
ylims!(1e-3,1e1)
xlims!(-6.0,-2.0)
# Sanity check seems to work! (At least at the level it did before for early time neutrinos)

k*bg.η[end]/Mpcfac
k*bg.η[end]/Mpcfac /36
η2x(bg.η[end] /12)

plot(bg.η(xx_k)*Mpcfac,𝒳₀_k[1](xx_k))
plot!(bg.η(xx_k)*Mpcfac,ie_0.s𝒳₀[1](xx_k),yscale=:log10)
plot!(perturb_conf.t,results_conf[2(ℓᵧ+1)+1,:])

hline!([0.0],color=:black)
vline!([3750])
exp(-η2x(3750))-1.0

plot!(bg.η(xx_k)*Mpcfac,𝒳₂_k[1](xx_k))
plot!(bg.η(xx_k)*Mpcfac,ie_0.s𝒳₂[1](xx_k))
plot!(perturb_conf.t,results_conf[2(ℓᵧ+1)+3,:])


results_conf
results
plot(bg.η(xx_k)*Mpcfac,𝒳₀_k[end](xx_k))
plot!(bg.η(xx_k)*Mpcfac,ie_0.s𝒳₀[end](xx_k))
plot!(perturb_conf.t,results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+n_q,:])
xlims!(0,2500)
ylims!(5,50)

#For this one we overshoot at the end if we use Ngamma3=500 or less, 700 is pretty good, more is diminishing returns
#Neutrinos have no effect
plot!(bg.η(xx_k)*Mpcfac,Θ₂(xx_k))
plot!(bg.η(xx_k)*Mpcfac,ie_0.sΘ2(xx_k))
plot!(perturb_conf.t,results_conf[3,:])
xlims!(0,2500)
ylims!(-.1,.1)


