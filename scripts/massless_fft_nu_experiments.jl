# Here we will start with the full hierarchy and switch to ie around horizon entry
using Bolt
using FFTW
using Plots
using DelimitedFiles
using Printf
using BenchmarkTools
using Interpolations
using Bolt: spline #FIXME why do I have to import this here but NOT in bg?
using LaTeXStrings
using OrdinaryDiffEq #TODO remove this when putting these functions into ie


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
# Background etc.
# 𝕡 = CosmoParams(
#     h = 0.6774,  # hubble factor
#     Ω_b = 0.0486, 
#     Ω_m = 0.2589,
#     Σm_ν = 0.00
# ) # Planck15 modifications to h, Ω_b,Ω_c, make mnu=0
𝕡 = CosmoParams(); 
n_q=15
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q);
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b); #FIXME γΩ
ih = IonizationHistory(𝕣, 𝕡, bg);
Mpcfac = bg.H₀*299792.458/100.
k = Mpcfac*kclass #get k in our units

# /// IC Free streaming ///
# Relevant Bessel functions (ℓ=0,1,2)
#ℓ=0
j0(x) = (x > 0.01) ? sin(x)/x : 1 - x^2 /6 + x^4 /120 - x^6 /5040
j0′(x) = -j1(x)
#ℓ=1
j1(x) =  (x > 0.01) ?  (sin(x) - x*cos(x))/x^2 : x/3 - x^3 /30 + x^5 /840
R1(x) =  (x > 0.01) ? j1(x) - 3j2(x)/x : 2x/15 - 2x^3 /105 + x^5 /1260
#ℓ=2
j2(x) = (x > 0.01) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^3 : x^2 /15 - x^4 /210 + x^6 /7560
j2′(x) = (x > 0.01) ? ( -x*(x^2 -9)*cos(x) + (4x^2 -9)*sin(x) ) / x^4 : 2x /15 - 2x^3 /105 + x^5 /1260
j2′′(x) = (x > 0.2) ? ( x*(5x^2 -36)*cos(x) + (x^4 - 17x^2 +36)*sin(x) ) / x^5 : 2/15 - 2x^2 /35 + x^4 /252 - x^6 /8910
R2(x) = (x > 0.2) ? -( j2(x) + 3j2′′(x) ) / 2 : -1/5 + 11x^2 /210 -x^4 /280 +17x^4 /166320
# The W coupling kernel (sum truncated at ℓ=2)
W00(x) = j0(x)
W01(x) = j1(x)
W02(x) = j2(x)
W21(x) = -R1(x)
W22(x) = -R2(x)
function Wsum(x,𝒳ᵢ₀,𝒳ᵢ₁,𝒳ᵢ₂)
    𝒳ₛ₀ = W00(x)*𝒳ᵢ₀ - 3W01(x)*𝒳ᵢ₁ + 5W02(x)*𝒳ᵢ₂  #ℓ=0 ( use the subscript ₛ for streaming, this is the "free-streaming" piece)
    𝒳ₛ₂ = W02(x)*𝒳ᵢ₀ - 3W21(x)*𝒳ᵢ₁ + 5W22(x)*𝒳ᵢ₂ #ℓ=2
    return 𝒳ₛ₀, 𝒳ₛ₂
end

# Hierarchy for comparison purposes - now replace with conformal hierarchy...
ℓᵧ=50
ℓ_mν=20
ℓ_ν=50#3#3#ℓ_ν10#ℓᵧ
pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
reltol=1e-12 
#solve the hierarchy just to be sure
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
results=zeros(pertlen,length(bg.x_grid))
perturb = boltsolve(hierarchy; reltol=reltol);
for (i_x, x) in enumerate(bg.x_grid)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
end

#conformal hierarchy
η2x = linear_interpolation(bg.η,bg.x_grid);
hierarchy_conf = ConformalHierarchy(hierarchy,η2x);
results_conf = boltsolve_conformal(hierarchy_conf;reltol=reltol);

plot(η2x.(results_conf.t/Mpcfac),results_conf(results_conf.t)[ν_idx,:])
plot!(bg.x_grid,results[ν_idx,:])
#WHAT IS HAPPENING HERE??

#truncated conformal hierarchy
# Input to the ie integrator struct (akin to hierarchy)
𝒩₀_0,𝒩₂_0 =  results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:] #hierarchy answer
spl0h𝒩₀,spl0h𝒩₂ = linear_interpolation(bg.x_grid,𝒩₀_0), linear_interpolation(bg.x_grid,𝒩₂_0)
ν_idx = 2(ℓᵧ+1) + 1
ie_0 = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
        spl0h𝒩₀,
        spl0h𝒩₂,
        ℓᵧ, ℓ_mν, n_q);
perturb_0 = boltsolve(ie_0;reltol=reltol); #no rsa

ie_0_conf = ConformalIEν(ie_0,η2x);
results_conf_ie_0 = boltsolve_conformal(ie_0_conf;reltol=reltol);

c𝒩₀_0,c𝒩₂_0 =  results_conf[2(ℓᵧ+1)+1,:],results_conf[2(ℓᵧ+1)+3,:] #hierarchy answer
c_spl0h𝒩₀,c_spl0h𝒩₂ = linear_interpolation(η2x(results_conf.t/Mpcfac),c𝒩₀_0), linear_interpolation(η2x(results_conf.t/Mpcfac),c𝒩₂_0)



#save  neutrinos:
# writedlm("./test/data/Bolt_mslss_nuperts_nonu_lmax$(ℓ_ν).dat",
#           hcat(bg.x_grid,results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:]))

#--- Begin neutrino functions ---#
function get_Φ′_Ψ(u,hierarchy::Hierarchy{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, par, bg, nq = hierarchy.k, hierarchy.par, hierarchy.bg,hierarchy.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩[2]
                                  + σℳ / bg.ρ_crit /4
                                  )
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩[0]
        + a^(-2) * ρℳ / bg.ρ_crit
        )
    return Φ′,Ψ
end

# Get the Φ' and Ψ (copy function in ie file) from hierarchy
function get_Φ′_Ψ(u,ie::IEν{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    𝒩₀ = ie.s𝒩₀(x)
    𝒩₂ = ie.s𝒩₂(x)
    ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2]+
                                  Ω_ν * 𝒩₂
                                  + σℳ / bg.ρ_crit /4
                                  )
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b
        + 4Ω_r * a^(-2) * Θ[0]
        + 4Ω_ν * a^(-2) * 𝒩₀
        + a^(-2) * ρℳ / bg.ρ_crit
        )
    return Φ′,Ψ
end

function fft_funcs(x, y, Φ′,Ψ, k,ℋ,q,m,𝕡)
    ϵ = (q^2 .+ exp.(2x)*m^2 ).^(1/2) #sqrt syntax doesn't work w/ bcast but put it back when undo bcast...
    q̃ = ϵ/q #convenience notation
    G₀ = ℋ .* q̃/k .* Φ′ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #for integrating in y #
    G₁ = -q̃.^2 .* Ψ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #
    K₀₀ = j0.(y) #1st index is ℓ 2nd index is derivative order
    K₀₁ = j0′.(y)
    K₂₀ = j2.(y) #
    K₂₁ = j2′.(y) #
    return G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁
end

function fft_integral(x, y,Φ′,Ψ,k,ℋ,q,m,𝕡,M)#,xd,yd,Φ′d,Ψd,ℋd) # for massive or massless neutrinos (𝒳=𝒩,ℳ)
    dy = y[2]-y[1]
    #  all ffts are performed in this function
    G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁ = fft_funcs(x,y, Φ′,Ψ, k,ℋ,q,m,𝕡) #
    # zero-pad the signals so convolution is not circular
    G₀,G₁ = [G₀; zeros(M-1)],[G₁; zeros(M-1)]
    K₀₀,K₀₁,K₂₀,K₂₁ = [K₀₀; zeros(M-1)],[K₀₁; zeros(M-1)],[K₂₀; zeros(M-1)],[K₂₁; zeros(M-1)] #
    # FFT the Gs, Ks
    G̃₀,G̃₁ = fft(G₀),fft(G₁)
    K̃₀₀, K̃₀₁, K̃₂₀, K̃₂₁ = fft(K₀₀),fft(K₀₁),fft(K₂₀),fft(K₂₁)#
    # Convolution theorem (iFFT pointwise product)
    𝒳₀ₓ = ifft(G̃₀.*K̃₀₀ .+ G̃₁.*K̃₀₁)[1:M]*dy 
    𝒳₂ₓ = ifft(G̃₀.*K̃₂₀ .+ G̃₁.*K̃₂₁)[1:M]*dy 
    return 𝒳₀ₓ,𝒳₂ₓ
end

function fft_ie(ie,perturb,M,m,q,i_q,u₀,x_grid)
    𝕡,bg,k,nq = ie_0.par,ie_0.bg,ie_0.k,ie.nq
    # Set up the "neutrino horizon" and FFT abscissas
    χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
    end
    𝒩₀ = ie.s𝒩₀(x_grid[1])
    𝒩₂ = ie.s𝒩₂(x_grid[1])
    _,_,𝒩, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
    if m==0 
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩,𝒩₂ )) #massless
    else
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
    end 
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
    # Put it all together
    # 𝒳₀,𝒳₂=zeros(Float32, M),zeros(Float32, M)
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)
end

function fft_ie_c(ie,perturb,M,m,q,i_q,u₀,x_grid)
    𝕡,bg,k,nq = ie_0.par,ie_0.bg,ie_0.k,ie.nq
    # Set up the "neutrino horizon" and FFT abscissas
    χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
    end
    𝒩₀ = ie.s𝒩₀(x_grid[1])
    𝒩₂ = ie.s𝒩₂(x_grid[1])
    _,_,𝒩, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
    if m==0 
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩,𝒩₂ )) #massless
    else
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
    end 
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
    # Put it all together
    # 𝒳₀,𝒳₂=zeros(Float32, M),zeros(Float32, M)
    # println("shapes, shape stream = $(length(𝒳ₛ₀)), shape fft = $(length(real.(𝒳₀ₓ)))")
    # println("M = $(M), $(length(yy)),$(length(invx))")
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)#,
end

function h_boltsolve_flex(hierarchy::Hierarchy{T},  x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T
    prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (x_ini , x_fin), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,
                dense=false,
                )
    return sol
end

function boltsolve_flex(ie::IEν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
    prob = ODEProblem{true}(Bolt.ie!, u₀, (x_ini , x_fin), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                # saveat=ie.bg.x_grid, 
                dense=false, #FIXME
                )
    return sol
end

function h_boltsolve_conformal_flex(confhierarchy::ConformalHierarchy{T},#FIXME we do't need this? {Hierarchy{T},AbstractInterpolation{T}},
    η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
    hierarchy = confhierarchy.hierarchy
    Mpcfac = hierarchy.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(Bolt.hierarchy_conformal!, u₀, 
                            (η_ini*Mpcfac , η_fin*Mpcfac),
                            confhierarchy)
    sol = solve(prob, ode_alg, reltol=reltol,
    dense=false
    )
    return sol
end

function boltsolve_conformal_flex(confie::ConformalIEν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
    η_ini,η_fin,u₀,ode_alg=KenCarp4(); reltol=1e-6) where T
    ie,η2x = confie.ie,confie.η2x
    Mpcfac = ie.bg.H₀*299792.458/100.
    prob = ODEProblem{true}(Bolt.ie_conformal!, u₀, 
                            (η_ini*Mpcfac, η_fin*Mpcfac),
                            confie)
    sol = solve(prob, ode_alg, reltol=reltol,
    dense=false
    )
    return sol
end


function iterate_fft(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
    M, reltol,x_ini, x_fin,u0) where T
    𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
    ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒩₀_km1, 𝒩₂_km1,
                    ℓᵧ, ℓ_mν, n_q)
    perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie(ie_k_late,perturb_k_late,M,0.,1.,0,
                        u0,perturb_k_late.t) #This is for massless only 
    return xx,𝒩₀_k,𝒩₂_k,perturb_k_late
end

function iterate_fft_c(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
    M, reltol,η_ini, η_fin,u0) where T
    𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
    ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒩₀_km1, 𝒩₂_km1,
                    ℓᵧ, ℓ_mν, n_q)
    ie_k_conf_late_c = ConformalIEν(ie_k_late,η2x);
    perturb_k_late_c = boltsolve_conformal_flex(ie_k_conf_late_c, η_ini, η_fin, u0; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie_c(ie_k_conf_late_c.ie,perturb_k_late_c,M,0.,1.,0,
                        u0,η2x(perturb_k_late_c.t/Mpcfac)) #This is for massless only 
    return xx,𝒩₀_k,𝒩₂_k,perturb_k_late_c
end

#---------------------------------#
# Itersolves
#---------------------------------#
function itersolve_fft(Nₖ::Int,ie_0::IEν{T},M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
    𝒩₀_0,𝒩₂_0 = ie_0.s𝒩₀,ie_0.s𝒩₂
    𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
    perturb_k = nothing
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft(𝒩₀_k,𝒩₂_k,ie_0.par,ie_0.bg,ie_0.ih,
                                   ie_0.k,ie_0.ℓ_γ,ie_0.ℓ_mν,ie_0.nq,
                                   M,reltol,x_ini,x_fin,u0)
    end
    return xx_k, 𝒩₀_k,𝒩₂_k,perturb_k
end
#ctime version
function itersolve_fft(Nₖ::Int,ie_0_c::ConformalIEν{T},M::Int,η_ini, η_fin,u0;reltol=1e-6) where T
    𝒩₀_0,𝒩₂_0 = ie_0_c.ie.s𝒩₀,ie_0_c.ie.s𝒩₂
    𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
    perturb_k = nothing
    ηη_k = nothing#zeros(length(bg.x_grid))
    for k in 1:Nₖ
        ηη_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft_c(𝒩₀_k,𝒩₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.ℓ_mν,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0)
    end
    return ηη_k,𝒩₀_k,𝒩₂_k,perturb_k
end



#---------------------------------#
#---------------------------------#
# Begin Experiments
#---------------------------------#
#---------------------------------#

function get_switch_u0(η,hierarchy_conf) #Input is η of the switch
    # switch_idx=740 #<- the switch idx for η=1.0ish
    hierarchy,bg = hierarchy_conf.hierarchy,hierarchy_conf.hierarchy.bg
    switch_idx = argmin(abs.(bg.η*Mpcfac .-η)) #for now we use the bg to find the switch
    #solve the split ode
    ℓᵧ = hierarchy.ℓᵧ
    # \/ we want to report this timing to get a full picture of total time (early+late)
    sol_early_c = h_boltsolve_conformal_flex(hierarchy_conf, bg.η[1], bg.η[switch_idx],  initial_conditions(bg.x_grid[1], hierarchy));
    # Get the new initial conditions
    u0_ie_c = zeros(2(ℓᵧ+1) + (0+1) + (ℓ_mν+1)*n_q + 5);
    for i in  1:2(ℓᵧ+1)
        u0_ie_c[i] = sol_early_c.u[end][i]
    end
    u0_ie_c[2(ℓᵧ+1)+1] = sol_early_c.u[end][2(ℓᵧ+1)+2]
    for i in  2(ℓᵧ+1)+(ℓ_ν+1)+1:pertlen
        down_shift = i-(ℓ_ν)
        u0_ie_c[down_shift] = sol_early_c.u[end][i]
    end
    return u0_ie_c
end


plot(bg.x_grid,c_spl0h𝒩₀.(bg.x_grid))
plot!(bg.x_grid,spl0h𝒩₀.(bg.x_grid))
plot!(bg.x_grid,results_conf(bg.η.(bg.x_grid)/Mpcfac)[ν_idx,:])
plot!(bg.x_grid,results[ν_idx,:])

# Set up the FFT struct an initial ansatz
# u0_ie_c[ν_idx]
# zero_ansatz₀,zero_ansatz₂ = linear_interpolation(η2x.(sol_late_c.t/Mpcfac),zeros(length(sol_late_c.t))), linear_interpolation(η2x.(sol_late_c.t/Mpcfac),zeros(length(sol_late_c.t)));
zero_ansatz₀,zero_ansatz₂ = linear_interpolation(η2x.(results_conf.t/Mpcfac),zeros(length(results_conf.t))), linear_interpolation(η2x.(results_conf.t/Mpcfac),zeros(length(results_conf.t)));
const_ansatz₀ = linear_interpolation(η2x.(results_conf.t/Mpcfac),u0_ie_c[ν_idx]*ones(length(results_conf.t)));
#^Shouldnn't matter
const_ansatz₀
zero_ansatz₀
# readdata = readdlm("./test/data/Bolt_mslss_nuperts_nonu_lmax3.dat")
readdata = readdlm("./test/data/Planck15_mslss_nuperts_nonu_lmax50.dat")
readx,read𝒩₀,read𝒩₂ = readdata[:,1],readdata[:,2],readdata[:,3]
ansatz₀,ansatz₂ = linear_interpolation(readx,read𝒩₀), linear_interpolation(readx,read𝒩₂);

ie_0_late_c = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    # zero_ansatz₀,
                    # const_ansatz₀,
                    # zero_ansatz₂,
                    # ansatz₀,
                    # ansatz₂,
                    spl0h𝒩₀,
                    spl0h𝒩₂,
                    ℓᵧ, ℓ_mν, n_q);
ie_0_conf_late_c = ConformalIEν(ie_0_late_c,η2x);


# Timings (@btime) to get the initial conditions via early full solve (k=0.3)
# 0.5 Mpc   182.678 ms (20831 allocations: 7.07 MiB)
# 1.0 Mpc   194.864 ms (21905 allocations: 7.10 MiB)
# 10.0 Mpc  263.569 ms (27421 allocations: 7.42 MiB)
# 100.0 Mpc 391.347 ms (36139 allocations: 8.13 MiB)

# And now ctime
M=2048*4
# Nᵢ=1
reltol=7e-4
#changing k, switch, hierarchy truncation, and ansatz will need to have a re-doing of u0_ie, ie_0 struct

# First experiment
η_switchη_switch = [10.0] #[0.5,1.0,10.0,100.0] 
η_switch = 10.0
MM = [8192]#[2^i for i in 12:14]
NᵢNᵢ = [2i-1 for i in 1:5] #max iters

#run this for plotting consistency
xx_kt,𝒩₀_kt,𝒩₂_kt,perturb_kt= itersolve_fft(1,ie_0_conf_late_c,MM[end],
    η_switchη_switch[1]/Mpcfac,bg.η[end],get_switch_u0(η_switchη_switch[1],hierarchy_conf);reltol=reltol);


# for η_switch in η_switchη_switch
    # Set the initial conditions at a particular switch value
    u0_ie_c = get_switch_u0(η_switch,hierarchy_conf)
    # Initial guess
    p1 = plot(xx_kt,ie_0_late_c.s𝒩₀.(xx_kt),label="I = 0, zero ansatz",legendfont=font(4),ls=:dash)
    plot!(p1,xx_kt,ie_0_late_c.s𝒩₀.(xx_kt),label="I = 0, mono init ansatz",legendfont=font(4),ls=:dash)
    plot!(p1,xx_kt,ie_0_late_c.s𝒩₀.(xx_kt),label="I = 0, Bolt init ansatz",legendfont=font(4),ls=:dash)
    plot!(p1,xx_kt,ie_0_late_c.s𝒩₀.(xx_kt),label="I = 0, Planck init ansatz",legendfont=font(4),ls=:dash)
    p2 = plot(xx_kt,ie_0_late_c.s𝒩₂.(xx_kt),label=false,c=p1.series_list[1][:linecolor],legendfont=font(4),ls=:dash)
    plot!(p2,xx_kt,ie_0_late_c.s𝒩₂.(xx_kt),label=false,legendfont=font(4),c=p1.series_list[end][:linecolor],ls=:dash)
    # Hierarchy
    plot!(p1,xx_kt,spl0h𝒩₀.(xx_kt),label="H",color=:black,lw=2)
    plot!(p2,xx_kt,spl0h𝒩₂.(xx_kt),label=false,color=:black,lw=2)
    plot!(p1,xx_kt,c_spl0h𝒩₀.(xx_kt),label="cspl",color=:red)
    plot!(p2,xx_kt,c_spl0h𝒩₂.(xx_kt),label=false,color=:red)


    for M in MM
        for Nᵢ in NᵢNᵢ
            @time xx_k,𝒩₀_k,𝒩₂_k,perturb_k = itersolve_fft(Nᵢ,ie_0_conf_late_c,M,
                                            η_switch/Mpcfac,bg.η[end],u0_ie_c;reltol=reltol);
            println("(M = $M, Nᵢ = $Nᵢ), 
                    error against full hierarchy is 
                    ℓ=0: $(sum( (spl0h𝒩₀.(xx_k) .- 𝒩₀_k.(xx_k)).^2 )/M),
                    ℓ=2: $(sum( (spl0h𝒩₂.(xx_k) .- 𝒩₂_k.(xx_k)).^2 )/M)\n")
            label="I = $(Nᵢ), Planck50_ansatz"#M = $(M), lmax = $(ℓ_ν)" 
            plot!(p1,xx_k,𝒩₀_k.(xx_k),label=label)#,c=:red)
            plot!(p2,xx_k,𝒩₂_k.(xx_k),label=false,c=p1.series_list[end][:linecolor])
                    
        end
    end
    # ylims!(p1,-0.1,0.1)
    # ylims!(p2,-0.05,0.04)
    xlabel!(p2,L"x",xguidefontsize=18)
    ylabel!(p1,L"\mathcal{N}_{0}",xguidefontsize=18)
    ylabel!(p2,L"\mathcal{N}_{2}",xguidefontsize=18)
    l = @layout [a  ; b]
    title!(p1,"k = $(@sprintf("%.2f", ie_0.k/Mpcfac
    )), vary ansatz, switch at $(@sprintf("%.1f", η_switch)) Mpc")
    p3 = plot(p1, p2, layout = l)
    savefig("../misc_plots/fft_debug/fft_experiments/mslss_k$(@sprintf("%.2f", ie_0.k/Mpcfac
            ))_switch$(@sprintf("%.1f", η_switch))_recheck_dipole_sc.pdf"
    )

# end

``
# x values of the horizon scale
η2x(3.4/Mpcfac)
η2x(37/Mpcfac)

η2x(10/Mpcfac)