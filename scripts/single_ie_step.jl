"""jms: Testing a single integral equation step"""
using Revise
using Bolt
using Plots
using NumericalIntegration
using Printf
using DelimitedFiles

# Load some saved hierarchy answers to compare against (and start from)
k_options = ["p03", "p3", "1p0"] #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
k_choice = k_options[1]
kMpc = parse(Float64, replace(k_choice,"p"=>"."))
ret = readdlm( @sprintf("./test/data/bolt_px_k%s.dat",k_choice) )
#for conveneince
dx = ret[2,1]-ret[1,1]
ics = ret[1,:]
step1 = ret[2,:]


# Write functions for the integrands of the photon IEs
#---

# Generalized version of our current g
function g(x̃, τ′,bg)
    τ_integrated = reverse(  cumul_integrate( reverse(x̃), reverse(τ′.(x̃)) )  )
    return @.(-τ′(x̃) * exp(-τ_integrated))
end

#check g
g(bg.x_grid[bg.x_grid.<bg.x_grid[end]] ,ih.τ′,bg)
plot(bg.x_grid[bg.x_grid.<bg.x_grid[end]],g(bg.x_grid[bg.x_grid.<bg.x_grid[end]],ih.τ′,bg))
plot!(bg.x_grid,ih.g̃(bg.x_grid)) #sanity check on g
xlims!(-7.5,-6)

#Temperature quadrupole integrand
function IΘ2(x, k,
              Π, Θ0, v_b, Φ′, Ψ,
              ih, bg)
    #all pert variables are at all x′ < x
    τ′,η = ih.τ′,bg.η #all splines of x
    x′= bg.x_grid[bg.x_grid.<x] #points do not include current timestep
    ḡ = g(x′ ,τ′, bg)
    y = @.(  k*( η(x)-η(x′) )  )#Bessel argument
    #NB we had to flip the signs here on all the collision terms to account for sign of our τ′
    IΘ2 = @.  ḡ*(  ( -Θ0 - Φ′/ τ′(x′) )*j2(y) + ( -v_b - k*Ψ / τ′(x′) )*j2′(y)  + Π*R2(y) / 2  )
    return IΘ2
end

#Polarization scalar integrand
function IΠ(x, k, Π, ih, bg)
    #Π is at all x′ < x
    x′= bg.x_grid[bg.x_grid.<x] #points do not include current timestep
    τ′,η = ih.τ′,bg.η #all splines of x
    ḡ = g(x′ ,τ′, bg)
    y = @.(  k*( η(x)-η(x′) )  )#Bessel argument
    IE2 = @. ḡ*j2bx2(y)*Π
    IΠ = 9IE2
    return IΠ
end

# Relevant Bessel functions #FIXME replace with Zack's code?
#TODO faster to spline these w/o losing sufficient accuracy?
#O(6) Taylor expansion expressions for small values of argument to avoid instability
#Here I am (almost) following Kamionkowski, but modify transtions slightly to reach minimum disagreement
#these are extremely accurate, 1e-14 tol in matching zone
#TODO move these to utils.jl
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
    #TODO: can streamline hierarchy and source funcs with this also
    k, ℓᵧ, par, bg, ih, nq = ie.k, 2, ie.par, ie.bg, ie.ih,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    ℓ_ν = ie.ℓ_ν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
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
#hierarchy for comparison purposes
function get_Φ′_Ψ(u, hierarchy::Hierarchy{T},x) where T
    k, ℓᵧ, par, bg, ih, nq = hierarchy.k, hierarchy.ℓᵧ, hierarchy.par, hierarchy.bg, hierarchy.ih,hierarchy.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    ℓ_ν = ie.ℓ_ν
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, hierarchy)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
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
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=2):ret[end,1], nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k = (bg.H₀*3e5/100)*kMpc

#test the ie integrator struct (akin to hierarchy)
ℓᵧ=2
ℓ_ν=50
ℓ_mν=20
reltol=1e-5 #cheaper  rtol
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(bg.x_grid))
ℳρ,ℳσ = zeros(length(bg.x_grid)),zeros(length(bg.x_grid)) #arrays for the massive neutrino integrated perts
ie = IE(BasicNewtonian(), 𝕡, bg, ih, k, 400, ℓ_ν, ℓ_mν, n_q)

#check ICs
ieic = initial_conditions(bg.x_grid[1],ie)
Φ′ic, Ψic = get_Φ′_Ψ(ieic,ie,bg.x_grid[1])
Θ0ic, v_bic = ieic[1],ieic[end]
Πic = Θ0ic + ieic[3] + ieic[4] + ieic[6]

#first temp step
IΘ2(bg.x_grid[1], k, Πic, Θ0ic, v_bic, Φ′ic, Ψic, ih, bg)
#This doesn't work, cumul_integrate in τ in g uses trapz, which wants at least 2 poitns
#^So we need at least a guess for what the value of Θ2 is after the first timestep
# or what all the other perturbation variables are after the first time step
#Can probably just use some analytic guess like we use to set up ICs? For one step early on not much happens anyways
#This first step is where Kamionkowski would use some TCA approximation (I think)

#punt on this and use hierarchy "answers" for now to check the integrand
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓ_ν, ℓ_ν, ℓ_mν,n_q)
all_extra_perts = hcat([ [p for p in get_Φ′_Ψ(ret[i,2:end],hierarchy,bg.x_grid[i])] for i in 1:length(ret[:,1]) ]...)
Φ′h,Ψh = all_extra_perts[1,:],all_extra_perts[2,:]

#check Ψ against Φ, and Φ' qualitatively
plot(bg.x_grid, log10.(-Φ′h))
plot(bg.x_grid, log10.(abs.(Ψh)))
plot!(bg.x_grid, log10.(abs.(ret[:,end-4])))

#plot integrated pert at each timestep using hierarchy "answer" for all previous steps
Θ0h = ret[:,1+1]
Πh, v_bh = ret[:,1+3] + ret[:,1+1+51]+ret[:,1+3+51], ret[:,end]
p0 = plot(ret[:,1],ret[:,1+3])
ieresT = [ integrate( ret[:,1][ret[:,1].<ret[i,1]],#bg.η(ret[:,1][ret[:,1].<ret[i,1]]).*bg.ℋ(ret[:,1][ret[:,1].<ret[i,1]]),
                     IΘ2(ret[i,1], k, Πh[1:i-1], Θ0h[1:i-1], v_bh[1:i-1], Φ′h[1:i-1], Ψh[1:i-1], ih, bg)
                    )
          for i in 4:length(ret[:,1])]
plot!(ret[4:end,1],ieresT )
xlims!(p0,-8,-3)
vline!(p0,[ret[1500,1]])
hline!(p0,[0])
#looks good!


#---
#Now to Π IE

#check pieces of Π
plot(ret[1:end,1],ret[1:end,1+3]+ret[1:end,1+1+51]+ret[1:end,1+3+51])
plot!(ret[1:end,1],ret[1:end,1+3])
plot!(ret[1:end,1],ret[1:end,1+2+50])
plot!(ret[1:end,1],ret[1:end,1+4+50]) #FIXME this looks bad...
xlims!(-8,-3)

#Check Π integrand in a regime where g is nonzero
p1 = plot()
istep=100
imin,imax=500,1700
for i in imin:istep:imax
    plot!(p1,ret[1:i,1],
         log10.(abs.(IΠ(ret[i,1], k, ret[1:i,1+3]+ret[1:i,1+1+51]+ret[1:i,1+3+51],ih,bg ))
         ))
end
p1

#the true evolution of Π from hierarchy - if equations are right should get the same thing back
Πh = ret[:,1+1+2]+ret[:,1+1+51]+ret[:,1+1+51+2]

#compare hierarchy Π against integrated Π in easy case of all previous steps being "correct"
plot(ret[:,1],Πh)
ieres = [ ret[i,1+3].+integrate( ret[:,1][ret[:,1].<ret[i,1]],#bg.η(ret[:,1][ret[:,1].<ret[i,1]]).*bg.ℋ(ret[:,1][ret[:,1].<ret[i,1]]),
                     IΠ(ret[i,1], k, Πh[1:i-1],ih,bg )
                    )
          for i in 4:length(ret[:,1])]#1500:1503]#length(ret[:,1]) ]

plot!(ret[4:end,1],ieres )
xlims!(-8,-3)
vline!([ret[1500,1]])
hline!([0])

#check the ratio
plot(ret[4:end,1],log10.(abs.(Πh[4:end] ./ ieres .-1)))
ret[:,1][ret[:,1].<=ret[1,1]]
#not great at early times, ok at middling/imporant times (but will want to improve), goes to zero at RSA times

#Π diagnostics - what is the bessel argument doing? the function? g?
p2,p3,p4=plot(),plot(),plot()
for i in 500:100:1700
    xtt = ret[:,1][i]
    x′t = ret[:,1][ret[:,1].<xtt]
    y = k.*( bg.η(xtt).-bg.η.(x′t) )
    plot!(p2,x′t,y)
    plot!(p3,x′t,j2bx2.(y))
    plot!(p4,x′t,j2bx2.(y).*  g(x′t ,ih.τ′, bg))
end
p2
xlims!(p2,-8,-3)

p3
xlims!(p3,-8,-3)

p4
xlims!(p4,-8,-3)

#---UNDER CONSTRUCTION

#More realistic test where "answers" are not provided

#Below is just pseudocode at the moment...

#Integrate the integrands (without the hierarchy crutch)
# function trapz(f,x) cumul_integrate()nd
#can  try cumul_integrate from NumericalIntegration.jl
xpts = ret[]
function θ2(x,x′ᵢ,N) cumul_integrate(Iθ2.(x′ᵢ:x:N,...),x′ᵢ:x:N) end
function Π(x, x′ᵢ,N) cumul_integrate(IΠ.(x′ᵢ:x:N,...),x′ᵢ:x:N) end
make sure the arguments are evaluated at all x' grid


#Picard iteration: ϕₖ₊₁(t) = u₀ + ∫₀ᵗdt′ f( t′, ϕₖ(t′) ), ϕ₀(t) = u₀
#FIXME would be better to do adaptive number of iters based on ϵ tol
function picard(u0,f,n_iter)
    for n in 1:n_iter
        ϕ₁ = u0 + θ2(x,x′ᵢ,N)
        ϕ₂ = u0 + Π(x,x′ᵢ,N)
    end

end

#matrix inversion?


# Try the iteration on IE

# Compare to hierarchy in a plot

# Test the iteration against saved hierarchy
# function test_iter()
#
#     @assert a == b
# end
