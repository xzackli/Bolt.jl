using Bolt
using FFTW
using Plots
using DelimitedFiles
using Printf
using NumericalIntegration
using Interpolations
using Bolt: spline #FIXME why do I have to import this here but NOT in bg?

# /// Setup ///

# Load some existing perturbations from the hierarchy
k_options = ["p03", "p3", "1p0"] #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
k_choice = k_options[1]
kMpc = parse(Float64, replace(k_choice,"p"=>".")) #DUMB does not work for comparison
ret = readdlm( @sprintf("./test/data/bolt_px_k%s_fine.dat",k_choice) )
retnf_class = open( @sprintf("./test/data/class_px_k%s_nofluid_re.dat",k_choice),"r" ) do datafile
# an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
#the second column is just a repeated k value, so remember it and delete col
kclass = retnf_class[2][1] #read class k mode from file (in h/Mpc)
dx = ret[2,1]-ret[1,1]
unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

# Background etc.
𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k = (bg.H₀*3e5/100)*kclass #get k in our units

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


# Hierarchy for comparison purposes
ℓ_ν=2
ℓᵧ=50
ℓ_mν=20
pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
reltol=1e-5 #cheaper  rtol
#solve the hierarchy just to be sure
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓᵧ, ℓ_mν,n_q)
results=zeros(pertlen+ℓᵧ-2,length(bg.x_grid))
perturb = boltsolve(hierarchy; reltol=reltol)
for (i_x, x) in enumerate(bg.x_grid)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
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
   
    𝒩[0] = ie.s𝒩₀(x)
    𝒩[2] = ie.s𝒩₂(x)#WHY DO WE NEED THIS HERE BUT NOT IN PHOTONS? AND NOT FOR MONO?

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

function fft_funcs(x, y, Φ′,Ψ, k,ℋ,q,m,𝕡)
    ϵ = (q^2 .+ exp.(2x)*m^2 ).^(1/2) #sqrt syntax doesn't work w/ bcast but put it back when undo bcast...
    q̃ = ϵ/q #convenience notation
    G₀ = ℋ .* q̃/k .* Φ′ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #for integrating in y
    G₁ = -q̃.^2 .* Ψ * (m==0. ? -1 : dlnf0dlnq(q,𝕡))
    K₀₀ = j0.(y) #1st index is ℓ 2nd index is derivative order
    K₀₁ = j0′.(y)
    K₂₀ = j2.(y)
    K₂₁ = j2′.(y)
    return G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁
end

function fft_integral(x, y,Φ′,Ψ,k,ℋ,q,m,𝕡) # for massive or massless neutrinos (𝒳=𝒩,ℳ)
    dy = y[2]-y[1]
    #  all ffts are performed in this function
    G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁ = fft_funcs(x,y, Φ′,Ψ, k,ℋ,q,m,𝕡)
    
    # zero-pad the signals so covolution is not circular
    M = length(y) 
    G₀,G₁ = [G₀; zeros(M-1)],[G₁; zeros(M-1)]
    K₀₀,K₀₁,K₂₀,K₂₁ = [K₀₀; zeros(M-1)],[K₀₁; zeros(M-1)],[K₂₀; zeros(M-1)],[K₂₁; zeros(M-1)]

    # FFT the Gs, Ks
    G̃₀,G̃₁ = fft(G₀),fft(G₁)
    K̃₀₀, K̃₀₁, K̃₂₀, K̃₂₁ = fft(K₀₀),fft(K₀₁),fft(K₂₀),fft(K₂₁)

    # Convolution theorem (iFFT pointwise product)
    𝒳₀ₓ = ifft(G̃₀.*K̃₀₀ .+ G̃₁.*K̃₀₁)[1:M]*dy
    𝒳₂ₓ = ifft(G̃₀.*K̃₂₀ .+ G̃₁.*K̃₂₁)[1:M]*dy
    return 𝒳₀ₓ,𝒳₂ₓ
end

function fft_ie(ie,perturb,M)
    x_grid = ie.bg.x_grid
    k = ie_0.k

    # Set up the "neutrino horizon" and FFT abscissas
    m0=0. #FIXME remove for massive neutrinos (in the future will just call the fft integral nq+1 times)
    q0=1. #It won't matter what q is since it will drop out of χν for massless, but set it to 1 
    χν0_   = spline([Bolt.χν(x, q0, m0 , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid], bg.x_grid)
    yy0x = k.*χν0_ 
    dy=(yy0x[end]-yy0x[1])/(M-1)
    yy0 = yy0x[1]:dy:yy0x[end]
    invx0 = LinearInterpolation(yy0x,x_grid).(yy0) #FIXME use spline?

    # Do the IC propagation
    u₀ = initial_conditions(first(x_grid), ie)
    _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)    
    𝒩ₛ₀, 𝒩ₛ₂ = unzip(Wsum.(yy0,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #test on massless

    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx0[j]),ie,invx0[j])
    end
   
    # Compute the new perts via FFT
    𝒩₀ₓ,𝒩₂ₓ = fft_integral(invx0, yy0, Φ′,Ψ, k, ie.bg.ℋ(invx0), q0,m0,ie.par) #massless

    # Put it all together
    𝒩₀ = 𝒩ₛ₀ .+ real.(𝒩₀ₓ) 
    𝒩₂ = 𝒩ₛ₂ .+ real.(𝒩₂ₓ)
    return invx0, LinearInterpolation(invx0,𝒩₀), LinearInterpolation(invx0,𝒩₂)

end

# Input to the ie integrator struct (akin to hierarchy)
𝒩₀_0,𝒩₂_0 =  results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:] #hierarchy answer
# 𝒩₀_0,𝒩₂_0 =  zeros(length(bg.x_grid)), zeros(length(bg.x_grid)) #correct answer
spl0h𝒩₀,spl0h𝒩₂ = spline(𝒩₀_0,bg.x_grid), spline(𝒩₂_0,bg.x_grid)
ie_0 = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
        spl0h𝒩₀,
        spl0h𝒩₂,
        length(bg.x_grid), ℓᵧ, ℓ_mν, n_q)
perturb_0 = boltsolve(ie_0;reltol=reltol) #no rsa

# run
M = 16*2048
xx0,res₀,res₂ = fft_ie(ie_0,perturb_0, M)
plot(xx0, res₀(xx0), label="ie-test-2048") #monopole
plot!(bg.x_grid,𝒩₀_0,label="hierarchy")

plot(ix0,res₂(ix0)) #quadrupole
plot!(bg.x_grid,𝒩₂_0)

println("test done")

#-------------------------------------------------------------------------------
#Debugging iteration below:

#FIXME rename these since calling it "integral equation" is confusing and wrong
# Do the iteration (this will be almost the same for massive)
function iterate(𝒩₀_km1,𝒩₂_km1, bg, ih, k, ℓᵧ, ℓ_mν, n_q,M)
    𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1)
    ie_k = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
            spline(𝒩₀_km1(bg.x_grid), bg.x_grid), #not ideal...
            spline(𝒩₂_km1(bg.x_grid), bg.x_grid),
            M, #FIXME this does nothing
            ℓᵧ, ℓ_mν, n_q)
    perturb_k = boltsolve(ie_k; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie(ie_k,perturb_k,M) 
    return xx,𝒩₀_k,𝒩₂_k
end

spline(𝒩₀_km1, xx)

xx1,res₀1,res₂1 = iterate(res₀,res₂,bg,ih,ie_0.k,ie_0.ℓ_γ,ie_0.ℓ_mν,ie_0.nq,M)

plot!(xx1, res₀1(xx1),label="ie-test-iter1")

#try manually...
𝒩₀_1,𝒩₂_1 = zero(res₀),zero(res₂)
ie_1 = IEν(BasicNewtonian(), 𝕡, ie_0.bg, ie_0.ih, ie_0.k,
    spline(res₀(bg.x_grid), bg.x_grid), #not ideal...
    spline(res₂(bg.x_grid), bg.x_grid),
    M, #FIXME this does nothing
    ie_0.ℓ_γ,ie_0.ℓ_mν,ie_0.nq)

ie_1

perturb_1 = boltsolve(ie_1; reltol=reltol)

#these splines don't match.
#Whe I pass the one-time FFT'd spline, I get garbage for the photon/neutrino monopole
plot(bg.x_grid, ie_1.s𝒩₀(bg.x_grid))
plot!(bg.x_grid, ie_0.s𝒩₀(bg.x_grid))
plot!(bg.x_grid,spline(res₀(bg.x_grid), bg.x_grid))
plot(bg.x_grid, ie_1.s𝒩₂(bg.x_grid))
plot!(bg.x_grid, ie_0.s𝒩₂(bg.x_grid))
plot!(bg.x_grid,spline(res₂(bg.x_grid), bg.x_grid))



plot(bg.x_grid,[perturb_0(bg.x_grid[i])[1] for i in 1:length(bg.x_grid)])
plot!(bg.x_grid,[perturb_1(bg.x_grid[i])[1] for i in 1:length(bg.x_grid)])

xx1,𝒩₀_1,𝒩₂_1 = fft_ie(ie_1,perturb_1,M) 
plot(𝒩₀_1.(xx1))

plot(bg.x_grid,[perturb_0(bg.x_grid[i])[2(ie_1.ℓ_γ+1)+1] for i in 1:length(bg.x_grid)])
plot!(bg.x_grid,[perturb_1(bg.x_grid[i])[2(ie_1.ℓ_γ+1)+1] for i in 1:length(bg.x_grid)])

plot(bg.x_grid,[perturb_0(bg.x_grid[i])[end-4] for i in 1:length(bg.x_grid)])
plot!(bg.x_grid,[perturb_1(bg.x_grid[i])[end-4] for i in 1:length(bg.x_grid)])
ylims!(0,2)

bg.η(-6.68) *(bg.H₀*299792.458/100  )
bg.η .* (bg.H₀*299792.458/100  )