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
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
results=zeros(pertlen,length(bg.x_grid));
perturb = boltsolve(hierarchy; reltol=reltol);
for (i_x, x) in enumerate(bg.x_grid)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
end
#conformal hierarchy
η2x = linear_interpolation(bg.η,bg.x_grid);

@time boltsolve(hierarchy; reltol=reltol);
@btime boltsolve(hierarchy; reltol=reltol);
#@btime 3.902 s (114143281 allocations: 2.47 GiB)
#FIXME do this for ctime

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

# Single iterations

# function iterate(Θ₂_km1,Π_km1, 𝕡::CosmoParams{T}, bg, ih, k, 
#     Nᵧ₁,Nᵧ₂,Nᵧ₃,xgi,
#     ℓ_ν, ℓ_mν, n_q,reltol) where T
#     Θ₂_k,Π_k = zero(Θ₂_km1),zero(Π_km1) #FIXME pre-allocate these (and below)
#     ie_k = IE(BasicNewtonian(), 𝕡, bg, ih, k,
#             linear_interpolation(xgi,Θ₂_km1),
#             linear_interpolation(xgi,Π_km1),
#             Nᵧ₁,Nᵧ₂,Nᵧ₃,
#             ℓ_ν, ℓ_mν, n_q)
#     u_all_k = boltsolve(ie_k; reltol=reltol)
#     N = length(xgi)
#     Φ′,Ψ,Θ₀,Π,v_b = zeros(N),zeros(N),zeros(N),zeros(N),zeros(N)
#     for (j,u) in enumerate( eachcol(u_all_k) )
#             Φ′[j],Ψ[j],Θ₀[j],Π[j],v_b[j] = get_perts(u,ie_k,xgi[j])
#     end
#     for i in 3:length(xgi)
#             Θ₂_k[i],Π_k[i] = g_weight_trapz_ie(i,xgi,ie_k,Φ′,Ψ,Θ₀,Π,v_b)
#     end
#     return Θ₂_k,Π_k,u_all_k
# end

# function iterate_fft_allν(𝒳₀_km1,𝒳₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, n_q,
#     M, reltol,x_ini, x_fin,u0) where T
#     ie_k_late = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
#                      𝒳₀_km1,𝒳₂_km1,
#                     ℓᵧ, n_q)
#     perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
#     xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k_late,perturb_k_late,M,
#                         u0,perturb_k_late.t) 
#     return xx,𝒳₀_k,𝒳₂_k,perturb_k_late
# end


#merged iterate
function iterate(Θ₂_km1,Π_km1, 𝒳₀_km1,𝒳₂_km1, 
                 ie::IEγν{T},
                #  𝕡::CosmoParams{T}, bg, ih, k, n_q,
                 M,x_ini, x_fin,u0,reltol) where T
    𝕡, bg, ih, k, n_q,Nᵧ₁,Nᵧ₂,Nᵧ₃ = ie.par,ie.bg,ie.ih,ie.k,ie.nq,ie.Nᵧ₁,ie.Nᵧ₂,ie.Nᵧ₃ #FIXME get rid of this line
    # println("photon types: ", typeof(Θ₂_km1),typeof(Π_km1))
    # println("neutrino types: ", typeof(𝒳₀_km1),typeof(𝒳₂_km1))
    # println("2: ", 𝒳₀_km1[1](-20.0))
    ie_k = IEγν(BasicNewtonian(), 𝕡, bg, ih, k,
                     Θ₂_km1,Π_km1,
                     𝒳₀_km1,𝒳₂_km1,
                     Nᵧ₁,Nᵧ₂,Nᵧ₃,
                     n_q)
    # println("3: ", 𝒳₀_km1[1](-20.0))
    # println("4: ", ie_k.s𝒳₀[1](-20.0))


    # Do the truncated boltzmann solve
    perturb_k = boltsolve_flex(ie_k, x_ini, x_fin, u0; reltol=reltol)
    
    # Get metric and photon-relevant perturbation variables
    #FIXME may want to put this in its own function
    xgi = x_grid_ie(ie_k,x_ini,x_fin)
    # println("xgi[1]",xgi[1])
    # println("x_ini,x_fin",x_ini,x_fin)
    N = length(xgi)
    # assert(length(xgi)==length(Π_k))
    Φ′,Ψ,Θ₀,Π,v_b = zeros(N),zeros(N),zeros(N),zeros(N),zeros(N)
    # for (j,u) in enumerate( eachcol(u_all_k) )
    #why is perturb_k allowed to be longer than N? probably because I don't enforce save_at_grid...
    # println("ptk size: ",size(perturb_k))
    u_all = perturb_k(xgi)
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
    xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k,M,u0,perturb_k.t,sΦ′,sΨ) 

    return xx,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k
end

# Merged itersolve
function itersolve(Nₖ::Int,ie_0::IEγν{T},M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
    𝒳₀_k,𝒳₂_k = ie_0.s𝒳₀,ie_0.s𝒳₂
    # println("1: ", 𝒳₀_k[1](-20.0))
    Θ₂_k,Π_k = ie_0.sΘ2,ie_0.sΠ
    xx_k,perturb_k = nothing,nothing
    for k in 1:Nₖ
        # println("iter = ",k)
        xx_k,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate(Θ₂_k,Π_k, 𝒳₀_k,𝒳₂_k,
                                                    ie_0,
                                                    M,x_ini,x_fin,u0,
                                                    reltol)
    end
    return xx_k,Θ₂_k,Π_k,𝒳₀_k,𝒳₂_k,perturb_k
end


# Helper functon for switch
function get_switch_u0(η,hierarchy) #Input is η of the switch
    # This function assumes truncated hierarchies for all neutrinos (but not yet photons)
    bg =hierarchy.bg
    Mpcfac = bg.H₀*299792.458/100.
    switch_idx = argmin(abs.(bg.η*Mpcfac .-η)) #for now we use the bg to find the switch
    #solve the split ode
    ℓᵧ,ℓ_ν,n_q = hierarchy.ℓᵧ,hierarchy.ℓ_ν, hierarchy.nq
    pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
    # \/ we want to report this timing to get a full picture of total time (early+late)
    sol_early_c = h_boltsolve_flex(hierarchy, bg.x_grid[1], bg.x_grid[switch_idx],  initial_conditions(bg.x_grid[1], hierarchy));
    
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
#^Is this really right? No cosmological modes of interest (up to k of 10) 
#enter the horizon before x=-13.3? which is 2π/k/300, ctime ~0.7 Mpc?
#^So for any reasonable k we should start before "horizon entry" so defined


# Initial setup
N_iters = 5 # in principle should replace this with a tolerance criterion
#sanity check (full hierarchy ansatz)
x_Θ₂_interp,x_Π₂_interp =  linear_interpolation(collect(bg.x_grid),results[3,:]),linear_interpolation(collect(bg.x_grid),results[3,:]+results[(ℓᵧ+1)+1,:]+results[(ℓᵧ+1)+3,:])
x_all_interps₀ = [linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+1,:]),[linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]) for idx_q in 1:n_q]...];
x_all_interps₂ = [linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+3,:]),[linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+idx_q,:]) for idx_q in 1:n_q]...];

ie_0 = IEγν(BasicNewtonian(), 𝕡, bg, ih, k,
        x_Θ₂_interp,x_Π₂_interp,
        x_all_interps₀,x_all_interps₂,
        50,50,400, #the first number is more imporant than I anticipated...
        n_q);
# constructor at least works initially, good...

η_switch = 1.0;
u0_ie = get_switch_u0(η_switch,hierarchy);
M = 2048*4

using Bolt

plot(bg.x_grid,x_all_interps₀[1].(bg.x_grid))

xx_k,Θ₂,Π,𝒳₀_k,𝒳₂_k,perturb_k = itersolve(N_iters,ie_0,M,η2x(η_switch/Mpcfac),0.0,u0_ie;reltol=reltol);


@btime itersolve(N_iters,ie_0,M,η2x(η_switch/Mpcfac),0.0,u0_ie;reltol=reltol);

plot(xx_k,Θ₂(xx_k))
plot!(xx_k,ie_0.sΘ2(xx_k))
plot(xx_k,Π(xx_k))
plot!(xx_k,ie_0.sΠ(xx_k))
plot(xx_k,𝒳₀_k[1](xx_k))
plot!(xx_k,ie_0.s𝒳₀[1](xx_k))
plot(xx_k,𝒳₂_k[1](xx_k))
plot!(xx_k,ie_0.s𝒳₂[1](xx_k))
plot(xx_k,abs.(𝒳₀_k[2](xx_k)),yscale=:log10)
plot!(xx_k,abs.(ie_0.s𝒳₀[2](xx_k)))
plot(xx_k,abs.(𝒳₂_k[2](xx_k)),yscale=:log10)
plot!(xx_k,abs.(ie_0.s𝒳₂[2](xx_k)))

# Sanity check seems to work! (At least at the level it did before for early time neutrinos)
