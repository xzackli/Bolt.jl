"""jms: Testing integral equation solver for single k mode"""
# using Revise - this doesn't work when running in vscode - revise I think is imported behind the scenes by default
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

# Relevant Bessel functions
#TODO faster to spline these w/o losing sufficient accuracy?
#O(6) Taylor expansion expressions for small values of argument to avoid instability
#Here I am (almost) following Kamionkowski, but modify transtions slightly to reach minimum disagreement
#these are extremely accurate, 1e-14 tol in matching zone
#TODO move these to utils.jl
# FIXME: Is ternary good or bad for the compiler?
j2(x) = (x > 0.01) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^3 : x^2 /15 - x^4 /210 + x^6 /7560
j2bx2(x)  =  (x > 0.06) ? -( 3x*cos(x) + (x^2 - 3)*sin(x) ) / x^5 : 1/15 - x^2 /210 + x^4 /7560 - x^6 /498960
j2′(x) = (x > 0.01) ? ( -x*(x^2 -9)*cos(x) + (4x^2 -9)*sin(x) ) / x^4 : 2x /15 - 2x^3 /105 + x^5 /1260
#these are less accurate, 1e-9ish tol in matching zone, have to match earlier
j2′′(x) = (x > 0.2) ? ( x*(5x^2 -36)*cos(x) + (x^4 - 17x^2 +36)*sin(x) ) / x^5 : 2/15 - 2x^2 /35 + x^4 /252 - x^6 /8910
R2(x) = (x > 0.2) ? -( j2(x) + 3j2′′(x) ) / 2 : -1/5 + 11x^2 /210 -x^4 /280 +17x^4 /166320



#function for generating metric inputs to IΘ2
#to save these inside the equivalent of hierarchy would confuse DE solver
#could hack and pass an ode with dy/dt = 0, but this is grosser than below
function get_Φ′_Ψ(u,ie::IE{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, ℓᵧ, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    ℓ_ν = ie.ℓ_ν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    # Θ[2] = ie.sΘ2(x)# call the spline, update Θ₂ at top since we do not evolve it
    ρℳ, σℳ  =  ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
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

# test the integrands
#---
#generate some background/ionization history
𝕡 = CosmoParams()
n_q=15
logqmin,logqmax = -6,-1
#different dx time resolutions - using slowest most accurate one to test
#^This is overkill and we don't need bg x_grid to be the same as ie x_grid
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q)
# bg = Background(𝕡; x_grid=-20.0:0.05:0.0, nq=n_q) # corase timesteps while testing
# bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=n_q) # corase timesteps while testing
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k = (bg.H₀*299792.458/100)*kclass #get k in our units

#input to the ie integrator struct (akin to hierarchy)
ℓᵧ=2
ℓ_ν=50
ℓ_mν=20
reltol=1e-8 #cheaper  rtol

#Kamionkowski weights
function Ws(xⱼ,xⱼ₊₁,τ,xᵢ)
    ϵτ = 1e-3 #if g is changing rapidly (τ′*dx>ϵτ), use g-aware weights #FIXME is the switch really necessary?
    dτ = -(τ(xⱼ₊₁) - τ(xⱼ))  #NB this is backwards from Kamionkowski since he does ``from 1''
    expτ = exp(-dτ)
    expj = exp( -( -τ(xᵢ) + τ(xⱼ₊₁) ) ) #NB ditto above
    τfac = (1 - (1+dτ)*expτ)/ dτ
    Wⱼ⁺ = expj* (  (dτ >ϵτ)  ? 1-expτ - τfac : dτ/2  )
    Wⱼ =  expj* (  (dτ >ϵτ)  ? τfac : dτ/2  )
    return Wⱼ,Wⱼ⁺
end

# KERNELS
function _IΘ2(x, x′,k,
    Π, Θ0, v_b, Φ′, Ψ,
    ih, bg) #for testing
    τ′,η = ih.τ′,bg.η #all splines of x
    y = @. k*( η(x)-η(x′) )#Bessel argument
    IΘ2 = @. ( ( Θ0 - Φ′/ (-τ′(x′))  )*j2(y) - 1 *( v_b   - ( k/bg.ℋ(x′) )*Ψ / (-τ′(x′)) )*j2′(y)  - 1 *Π*R2(y) / 2  )
    return IΘ2
end

function _IΠ(x, x′,k, Π, ih, bg)
    η = bg.η #all splines of x
    y = @. k*( η(x)-η(x′) )#Bessel argument
    IE2 = @. j2bx2(y)*Π
    IΠ = 9IE2
    return IΠ
end

#Without trying anything new...
function g_weight_trapz_ie(x,x_grid,ie,u_all)
    # x_grid = ie.bg.x_grid
    τ = ie.ih.τ
    k = ie.k
    Θ,Θᴾ,𝒩,ℳ,Φ,δ,v,δ_b,v_b = ie_unpack(u_all,ie)
    #update Θ₂, Π with splines at all x by mutation-> first set Θ₂, then use Π, hierarchy Θ₀ᴾ to get Θ₂ᴾ
    Θ[2,:] .= ie.sΘ2.(x_grid)
    Π = ie.sΠ.(x_grid)#@.Θ[2] + Θᴾ[0] + Θᴾ[2]
    Θᴾ[2,:] .= Π .- Θᴾ[0,:] .- Θ[2,:]
    #probably not optimal to do this here at every step?...
    Φ′,Ψ = zeros(length(x_grid)),zeros(length(x_grid))
    for (j,u) in enumerate( eachcol(u_all) )
        # println("j = ", j)
     Φ′[j],Ψ[j] = get_Φ′_Ψ(u,ie,x_grid[j])
    end
    #TODO can just dot notation this whole operation?
    #TODO or use cumul integrate? would it be faster?
    i = length(x_grid[x_grid.<=x])
    xᵢ = x
    Θ2ᵢ, Πᵢ = 0,0 
    Wⱼ,Wⱼ⁺ = 0, 0 
    for j in 1:i-2
        xⱼ,xⱼ₊₁ = x_grid[j], x_grid[j+1]
        Wⱼ,Wⱼ⁺ = Ws(xⱼ,xⱼ₊₁,τ,xᵢ) #passing xᵢ for now but could update later externally...
        #TODO if we want to compute weight once and save them we can?
        #Implicit weighting scheme at each timestep
        Θ2ᵢ += (_IΘ2(xᵢ,xⱼ₊₁,k,Π[j+1],Θ[0,j+1],v_b[j+1],Φ′[j+1],Ψ[j+1],ie.ih,ie.bg)*Wⱼ⁺
               + _IΘ2(xᵢ,xⱼ,k,Π[j],Θ[0,j],v_b[j],Φ′[j],Ψ[j],ie.ih,ie.bg)*Wⱼ)
        Πᵢ += ( _IΠ(xᵢ,xⱼ₊₁,k,Π[j+1],ie.ih,ie.bg)*Wⱼ⁺
               + _IΠ(xᵢ,xⱼ,k,Π[j],ie.ih,ie.bg)*Wⱼ)
    end
    #Handle final sub-timestep j = i-1 (pull out final loop iteration)
    #FIXME this is silly just put an if in the loop?
    xⱼ,xⱼ₊₁ = x_grid[i-1], x_grid[i]
    Wⱼ,Wⱼ⁺ = Ws(xⱼ,xⱼ₊₁,τ,xᵢ) #passing xᵢ for now but could update later externally...
    Θ2ᵢ += (_IΘ2(xᵢ,x_grid[i],k,0.,Θ[0,i],v_b[i],Φ′[i],Ψ[i],ih,bg)*Wⱼ⁺
           + _IΘ2(xᵢ,x_grid[i-1],k,Π[i-1],Θ[0,i-1],v_b[i-1],Φ′[i-1],Ψ[i-1],ih,bg)*Wⱼ)
    Πᵢ += _IΠ(xᵢ,x_grid[i-1],k,Π[i-1],ih,bg)*Wⱼ
    #Kamionkowski integration scheme for handling x′ = x at each x (the implicit timestep)
    #Not sure how exactly derived, but apparently comes from the matrix inversion
    Πᵢ = (Πᵢ + Θ2ᵢ) / ( 1 - 7/10 * Wⱼ⁺)
    Θ2ᵢ = Θ2ᵢ + ( i <length(x_grid) ? Π[i+1] : 0. )/10 * Wⱼ⁺ #if i+1>length(x_grid), return 0 for oob array
    return Θ2ᵢ,Πᵢ
end


#FIXME consolidate on interpolator -> pass an interpolator rather than the ingredients
function iterate(Θ₂_km1,Π_km1, bg, ih, k, 
                 Nᵧ₁,Nᵧ₂,Nᵧ₃,xgi,
                 ℓ_ν, ℓ_mν, n_q,reltol)
    Θ₂_k,Π_k = zero(Θ₂_km1),zero(Π_km1)
    ie_k = IE(BasicNewtonian(), 𝕡, bg, ih, k,
            linear_interpolation(xgi,Θ₂_km1),
            linear_interpolation(xgi,Π_km1),
            Nᵧ₁,Nᵧ₂,Nᵧ₃,
            ℓ_ν, ℓ_mν, n_q)
    u_all_k = boltsolve(ie_k; reltol=reltol)
#     print("lxgi: ",length(xgi), " U: ", length(u_all_k[:,1]))
    for i in 3:length(xgi)
        # println("i= ",i)
        Θ₂_k[i],Π_k[i] = g_weight_trapz_ie(xgi[i],xgi,
                                           ie_k,u_all_k) 
    end
    return Θ₂_k,Π_k,u_all_k
end
function iterate_spl(Θ₂_km1,Π_km1, bg, ih, k, 
        Nᵧ₁,Nᵧ₂,Nᵧ₃,xgi,
        ℓ_ν, ℓ_mν, n_q,reltol)
        Θ₂_k,Π_k = zero(Θ₂_km1),zero(Π_km1)
        ie_k = IE(BasicNewtonian(), 𝕡, bg, ih, k,
                        spline(Θ₂_km1, xgi),
                        spline(Π_km1, xgi),
                        Nᵧ₁,Nᵧ₂,Nᵧ₃,
                        ℓ_ν, ℓ_mν, n_q)
                
        u_all_k = boltsolve(ie_k; reltol=reltol)
        for i in 3:length(xgi)
                Θ₂_k[i],Π_k[i] = g_weight_trapz_ie(xgi[i],xgi,
                                                ie_k,u_all_k) 
        end
        return Θ₂_k,Π_k,u_all_k
end



# uncollect(v) = v[1]:(v[end]-v[1])/(length(v)-1):v[end]

#initialize splines to zero
xgrid_hier = ret[1,1]:round(ret[2,1]-ret[1,1],digits=3):ret[end,1]
Θ₂_0,Π_0 = zeros(length(xgrid_hier)),zeros(length(xgrid_hier))
# Θ₂_0,Π_0 = 1. *ret[:,4], 1. * ( 1. *ret[:,1+3] .+ 1. *(ret[:,1+1+51] .+ ret[:,1+3+51]) )
#^This is here to test what happens if you give the right answer on the first iter (nothing, which is what should happen)
# spl0hΘ₂,spl0hΠ = spline(Θ₂_0,xgrid_hier), spline(Π_0,xgrid_hier)
spl0hΘ₂,spl0hΠ = spline(Θ₂_0,xgrid_hier), spline(Π_0,xgrid_hier) #FIXME names 
#initialize an IE object for the first iter...
#FIXME: when you declare an IE struct the fields are immutable? So can't update splines...new object for each iter...probably bad
indhor = argmin(abs.(k .* bg.η .- 1))
indrec = argmin(abs.( -ih.τ′ .* bg.ℋ .*bg.η .- 1))
xhor,xdec = bg.x_grid[indhor],bg.x_grid[indrec]
ie_0 = IE(BasicNewtonian(), 𝕡, bg, ih, k,
        spl0hΘ₂,
        spl0hΠ,
        indhor,indrec-indhor,length(bg.x_grid)-indrec,
        # these values are for the test case, in reality we won't need them after testing
        # this doesnt work, was trying to make like original xgrid length(bg.x_grid)÷3, length(bg.x_grid)÷3, length(bg.x_grid)÷3, 
        ℓ_ν, ℓ_mν, n_q)
#do the first ODE solve
u_all_0 = boltsolve(ie_0; reltol=reltol)#_rsa(ie_0)


#try 1 iter
Θ₂_0,Π_0 = spl0hΘ₂.(ie_0.bg.x_grid), spl0hΠ.(ie_0.bg.x_grid)
Θ₂_1,Π_1,u_all_1 = iterate_spl(Θ₂_0,Π_0, ie_0.bg, ie_0.ih, ie_0.k,  
                      indhor,indrec-indhor,length(bg.x_grid)-indrec, bg.x_grid,ℓ_ν, ℓ_mν, n_q,
                      reltol)
println(indhor, ", ", indrec-indhor, ", ",length(bg.x_grid)-indrec)
# Try a sparser set of xgrid points

# We can leave the bg xgrid points as before, but now use the supplied x points in IE
# Need to modify the iter source for this...do NOT just get the bg xgrid
# We can either make ie take x_points as a struct field, or pass that to trapz_weights
# The easiest thing is the former, since then wherever ie goes the x_grid goes...
x_grid_ie(ie_0) ≈ collect(bg.x_grid) #check that we reproduce the old (slow/fine) grid



litp0hΘ₂,litp0hΠ = linear_interpolation(xgrid_hier,Θ₂_0), linear_interpolation(xgrid_hier,Π_0) #FIXME names 
# ie_0_sparse.Nᵧ₁,ie_0_sparse.Nᵧ₂,ie_0_sparse.Nᵧ₃
ie_0_sparse = IE(BasicNewtonian(), 𝕡, bg, ih, k,
        litp0hΘ₂,
        litp0hΠ,
        20,100,200,
        # these values are for the test case, in reality we won't need them after testing
        # this doesnt work, was trying to make like original xgrid length(bg.x_grid)÷3, length(bg.x_grid)÷3, length(bg.x_grid)÷3, 
        ℓ_ν, ℓ_mν, n_q)
#do the first ODE solve
u_all_0_sparse = boltsolve(ie_0_sparse; reltol=reltol)
#try 1 iter
x_grid_ie_0_sparse = x_grid_ie(ie_0_sparse)
Θ₂_0_sparse,Π_0_sparse = litp0hΘ₂.(x_grid_ie_0_sparse), litp0hΠ.(x_grid_ie_0_sparse)

Θ₂_1_sparse,Π_1_sparse,u_all_1_sparse = iterate(Θ₂_0_sparse,Π_0_sparse, ie_0_sparse.bg, ie_0_sparse.ih, ie_0_sparse.k,  
                                 ie_0_sparse.Nᵧ₁,ie_0_sparse.Nᵧ₂,ie_0_sparse.Nᵧ₃, 
                                 x_grid_ie_0_sparse, ℓ_ν, ℓ_mν, n_q, reltol)

#more sparse iters
Θ₂_2_sparse,Π_2_sparse,u_all_2_sparse = iterate(Θ₂_1_sparse,Π_1_sparse, ie_0_sparse.bg, ie_0_sparse.ih, ie_0_sparse.k,  
                                 ie_0_sparse.Nᵧ₁,ie_0_sparse.Nᵧ₂,ie_0_sparse.Nᵧ₃, 
                                 x_grid_ie(ie_0_sparse),
                                 ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_3_sparse,Π_3_sparse,u_all_3_sparse = iterate(Θ₂_2_sparse,Π_2_sparse, bg, ih, k, 
                                 ie_0_sparse.Nᵧ₁,ie_0_sparse.Nᵧ₂,ie_0_sparse.Nᵧ₃, 
                                 x_grid_ie(ie_0_sparse),
                                 ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_4_sparse,Π_4_sparse,u_all_4_sparse = iterate(Θ₂_3_sparse,Π_3_sparse, bg, ih, k, 
                                     ie_0_sparse.Nᵧ₁,ie_0_sparse.Nᵧ₂,ie_0_sparse.Nᵧ₃, 
                                     x_grid_ie(ie_0_sparse),
                                     ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_5_sparse,Π_5_sparse,u_all_5_sparse2 = iterate(Θ₂_4_sparse,Π_4_sparse, bg, ih, k, 
                                     ie_0_sparse.Nᵧ₁,ie_0_sparse.Nᵧ₂,ie_0_sparse.Nᵧ₃, 
                                     x_grid_ie(ie_0_sparse),
                                     ℓ_ν, ℓ_mν, n_q, reltol)




#preliminary timing

#more iters - (I have a function for this, but manually do it so we can check results)
Θ₂_0,Π_0 = spl0hΘ₂.(ie_0.bg.x_grid), spl0hΠ.(ie_0.bg.x_grid)
Θ₂_1,Π_1,u_all_1 = iterate_spl(Θ₂_0,Π_0, ie_0.bg, ie_0.ih, ie_0.k,  
                      indhor,indrec-indhor,length(bg.x_grid)-indrec, bg.x_grid,ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_2,Π_2,u_all_2 = iterate_spl(Θ₂_1,Π_1, bg, ih, k, 
                        ie_0.Nᵧ₁,ie_0.Nᵧ₂,ie_0.Nᵧ₃, bg.x_grid,
                        ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_3,Π_3,u_all_3 = iterate_spl(Θ₂_2,Π_2, bg, ih, k, 
                   ie_0.Nᵧ₁,ie_0.Nᵧ₂,ie_0.Nᵧ₃, bg.x_grid,
                   ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_4,Π_4,u_all_4 = iterate_spl(Θ₂_3,Π_3, bg, ih, k, 
                       ie_0.Nᵧ₁,ie_0.Nᵧ₂,ie_0.Nᵧ₃, bg.x_grid,
                       ℓ_ν, ℓ_mν, n_q, reltol)
Θ₂_5,Π_5,u_all_5 = iterate_spl(Θ₂_4,Π_4, bg, ih, k, 
                       ie_0.Nᵧ₁,ie_0.Nᵧ₂,ie_0.Nᵧ₃, bg.x_grid,
                       ℓ_ν, ℓ_mν, n_q, reltol)



#check hierarchy & timesteps
η2x = linear_interpolation(bg.η,bg.x_grid)
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, 50, ℓ_ν, ℓ_mν,n_q)
results = boltsolve(hierarchy;reltol=1e-5)

hierarchy_conf = ConformalHierarchy(hierarchy,η2x)
results_conf = boltsolve_conformal(hierarchy_conf;reltol=1e-5)


#temp mono
plot(ret[:,1],ret[:,1+1],label="hierarchy",color=:red)
plot!(results_conf.t,results_conf[1,:],label="hierarchy-conf",color=:orange,ls=:dash,xscale=:log10)
plot(bg.η.(results.t)*bg.H₀*3e5/100,results[1,:],label="hierarchy-updated-rt1e-5",color=:blue)
xlims!(bg.η(-12)*bg.H₀*3e5/100,bg.η(0)*bg.H₀*3e5/100)

plot!(bg.η.(retnf_class[1,:])*bg.H₀*3e5/100,-retnf_class[2+1,:]/4/𝕡.h,label="hierarchy-CLASS",color=:green,ls=:dash)
plot(bg.x_grid,u_all_0[1,:],label="iter 0",ls=:dot)
vline!([bg.x_grid[indrec]],ls=:dash,color=:black)

plot(bg.x_grid,u_all_1[1,:],ls=:dash,label="iter 1")
plot(bg.x_grid,u_all_2[1,:],ls=:dash,label="iter 2")
plot!(bg.x_grid,u_all_3[1,:],ls=:dash,label="iter 3")
plot!(bg.x_grid,u_all_4[1,:],ls=:dash,label="iter 4")
plot!(bg.x_grid,u_all_5[1,:],ls=:solid,label="iter 5 - global",color=:black)
plot!(x_grid_ie_0_sparse, u_all_0_sparse[1,:],label="sparse-iter 0",ls=:dot)
plot!(x_grid_ie_0_sparse, u_all_1_sparse[1,:],label="sparse-iter 1",ls=:dash)
plot!(x_grid_ie_0_sparse,u_all_2_sparse[1,:],ls=:dash,label="sparse-iter 2")
plot!(x_grid_ie_0_sparse,u_all_3_sparse[1,:],ls=:dash,label="sparse-iter 3")
plot!(x_grid_ie_0_sparse,u_all_4_sparse[1,:],ls=:dash,label="sparse-iter 4")
plot!(x_grid_ie_0_sparse,u_all_5_sparse[1,:],ls=:dash,label="sparse-iter 5")
xlims!(-12,0)
vline!([xhor,xdec],ls=:dash,label="transitions")
# ylims!(-.5,.5)
xlabel!("x")
ylabel!("temp mono")
savefig("../temp_mono_hier_N3-200-2000fix_rtol-8.png")

#temp quadrupole

plot(ret[:,1],ret[:,1+3],label="hierarchy",color=:red)
plot(x_grid_ie_0_sparse,Π ,ls=:dash,label="sparse-iter 0")
plot!(x_grid_ie_0_sparse,u_all_1_sparse[3,:],ls=:dash,label="sparse-iter 1")
plot!(x_grid_ie_0_sparse,u_all_2_sparse[3,:],ls=:dash,label="sparse-iter 2")
plot!(x_grid_ie_0_sparse,u_all_3_sparse[3,:],ls=:dash,label="sparse-iter 3")
plot!(x_grid_ie_0_sparse,u_all_4_sparse[3,:],ls=:dash,label="sparse-iter 4")
plot!(x_grid_ie_0_sparse,@.(Π_5_sparse - u_all_5_sparse[3+1,:]- u_all_5_sparse[3+3,:]),ls=:dash,label="sparse-iter 5")
plot!(bg.x_grid,@.(Π_5- u_all_5[3+1,:]- u_all_5[3+3,:]),ls=:dash,label="iter 5")




#polzn mono
plot(ret[:,1],ret[:,1+51+1],label="hierarchy",color=:black)
plot!(bg.x_grid,u_all_0[3+1,:],ls=:dot,label="iter 0")
plot!(bg.x_grid,u_all_1[3+1,:],ls=:dash,label="iter 1")
plot!(bg.x_grid,u_all_2[3+1,:],ls=:dash,label="iter 2")
plot!(bg.x_grid,u_all_3[3+1,:],ls=:dash,label="iter 3")
plot!(bg.x_grid,u_all_4[3+1,:],ls=:dash,label="iter 4")
plot!(bg.x_grid,u_all_5[3+1,:],ls=:dash,label="iter 5")
ylims!(-0.005,0.017)
xlims!(-9,0)

#polzn dipole
plot(ret[:,1],ret[:,1+51+2],label="hierarchy",color=:black)
plot!(bg.x_grid,u_all_0[3+2,:],ls=:dot,label="iter 0")
plot!(bg.x_grid,u_all_1[3+2,:],ls=:dash,label="iter 1")
plot!(bg.x_grid,u_all_2[3+2,:],ls=:dash,label="iter 2")
plot!(bg.x_grid,u_all_3[3+2,:],ls=:dash,label="iter 3")
plot!(bg.x_grid,u_all_4[3+2,:],ls=:dash,label="iter 4")
plot!(bg.x_grid,u_all_5[3+2,:],ls=:dash,label="iter 5")
ylims!(-0.005,0.017)
xlims!(-9,0)


