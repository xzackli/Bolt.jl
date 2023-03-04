# Here we will start with the full hierarchy and switch to ie around horizon entry
using Bolt
using FFTW
using Plots
using Plots.PlotMeasures
using DelimitedFiles
using Printf
using BenchmarkTools
using Interpolations
using Bolt: spline #FIXME why do I have to import this here but NOT in bg?
using LaTeXStrings
using OrdinaryDiffEq #TODO remove this when putting these functions into ie


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
function get_Φ′_Ψ(u,ie::IEallν{T},x) where T
    #TODO: can streamline hierarchy and source funcs with this helper function also
    k, par, bg, nq = ie.k, ie.par, ie.bg,ie.nq
    Ω_r, Ω_b, Ω_m, N_ν, H₀² = par.Ω_r, par.Ω_b, par.Ω_m, par.N_ν, bg.H₀^2 #add N_ν≡N_eff
    ℋₓ =  bg.ℋ(x)
    a = x2a(x)
    Ω_ν =  7*(2/3)*N_ν/8 *(4/11)^(4/3) *Ω_r
    Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(u, ie)  # the Θ, Θᵖ, 𝒩 are views (see unpack)
    𝒩[0] = ie.s𝒳₀[1](x)
    𝒩[2] = ie.s𝒳₂[1](x)
    for idx_q in 0:(nq-1)
        ℳ[0*nq+idx_q] = ie.s𝒳₀[idx_q+1](x)
        ℳ[2*nq+idx_q] = ie.s𝒳₂[idx_q+1](x)
    end
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
    G₀ = ℋ .* q̃/k .* Φ′ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #for integrating in y #
    G₁ = -q̃.^2 .* Ψ * (m==0. ? -1 : dlnf0dlnq(q,𝕡)) #
    K₀₀ = j0.(y) #1st index is ℓ 2nd index is derivative order
    K₀₁ = j0′.(y)
    K₂₀ = j2.(y) #
    K₂₁ = j2′.(y) #
    return G₀,K₀₀,K₀₁, G₁,K₂₀,K₂₁
end

function fft_integral(x, y,Φ′,Ψ,k,ℋ,q,m,𝕡,M) # for massive or massless neutrinos (𝒳=𝒩,ℳ)
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

function fft_ie(ie::IEν,perturb,M,m,q,i_q,u₀,x_grid)
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    # Set up the "neutrino horizon" and FFT abscissas
    # χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    χνs = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q1,𝕡.Σm_ν) for x in bg.x_grid])
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
    end
    _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
    if m==0 
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
    else
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
    end 
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)
end

function fft_ie_c(ie::IEν,perturb,M,m,q,i_q,u₀,x_grid) #FIXME add type decorators
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    # Set up the "neutrino horizon" and FFT abscissas
    # χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    χνs = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q1,𝕡.Σm_ν) for x in bg.x_grid])
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
    end
    _,_,𝒩₀, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
    if m==0 
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
    else
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
    end 
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,m,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    return invx, linear_interpolation(invx,𝒳₀), linear_interpolation(invx,𝒳₂)#,
end


function fft_ie(ie::IEallν,perturb,M,u₀,x_grid)
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    # Set up the "neutrino horizon" and FFT abscissas
    # χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    #empty splines
    all_splines₀ = Array{AbstractInterpolation}(undef,nq+1)
    all_splines₂ = Array{AbstractInterpolation}(undef,nq+1)

    #explicitly do massless case
    χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),1.0,0.0,bg.quad_pts,bg.quad_wts) for x in x_grid]) #bg.η
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources #FIXME this should probably happen outside of this function
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
    end
    _,_,𝒩₀,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), 1.0,0.0,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    all_splines₀[1] = linear_interpolation(invx,𝒳₀)
    all_splines₂[1] = linear_interpolation(invx,𝒳₂)

    #massive case
    for i_q in 0:nq-1
        q = Bolt.xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),q,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts) for x in x_grid])
        yyx = k.* (χνs .- χνs[1])
        dy=(yyx[end]-yyx[1])/(M-1)
        yy = yyx[1]:dy:yyx[end]
        invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
        
        # Get metric sources #FIXME this should probably happen outside of this function
        Φ′,Ψ = zeros(M),zeros(M)
        for j in 1:M
            Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
        end
        _,_,_, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
        
        # Compute the new perts via FFT
        𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,Σm_ν,𝕡,M)#,
        # Put it all together
        𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
        𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 

        all_splines₀[i_q+2] = linear_interpolation(invx,𝒳₀)
        all_splines₂[i_q+2] = linear_interpolation(invx,𝒳₂)

    end


    return invx, all_splines₀, all_splines₂
end

function fft_ie_c(ie::IEallν,perturb,M,u₀,x_grid) #FIXME add type decorators
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    # Set up the "neutrino horizon" and FFT abscissas
    # χνs = [Bolt.χν(x, q, m , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in x_grid]
    #empty splines
    # all_splines₀ = [nothing for i in 1:nq+1] # doesnlt work
    all_splines₀ = Array{AbstractInterpolation}(undef,nq+1)
    all_splines₂ = Array{AbstractInterpolation}(undef,nq+1)

    #explicitly do massless case
    χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),1.0,0.0,bg.quad_pts,bg.quad_wts) for x in x_grid])
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    println("χνs : ", χνs)
    println("dy: ", dy)
    yy = yyx[1]:dy:yyx[end]
    println("yy shape: ", size(yy))
    println("yyx shape: ", size(yyx))
    println("χνs shape: ", size(χνs))
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources #FIXME this should probably happen outside of this function
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
    end
    _,_,𝒩₀,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀[0],𝒩₀[1],𝒩₀[2])) #massless
    # Compute the new perts via FFT
    𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), 1.0,0.0,𝕡,M)#,
    # Put it all together
    𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
    𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 
    all_splines₀[1] = linear_interpolation(invx,𝒳₀)
    all_splines₂[1] = linear_interpolation(invx,𝒳₂)

    #massive case
    for i_q in 0:nq-1
        q = Bolt.xq2q(bg.quad_pts[i_q+1],logqmin,logqmax)
        χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),q,𝕡.Σm_ν,bg.quad_pts,bg.quad_wts) for x in x_grid])
        yyx = k.* (χνs .- χνs[1])
        dy=(yyx[end]-yyx[1])/(M-1)
        yy = yyx[1]:dy:yyx[end]
        invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
        
        # Get metric sources #FIXME this should probably happen outside of this function
        Φ′,Ψ = zeros(M),zeros(M)
        for j in 1:M
            Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb(invx[j]),ie,invx[j])
        end
        _,_,_, ℳ₀,_,_,_,_,_ =  unpack(u₀,ie)   
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q])) #massive
        
        # Compute the new perts via FFT
        𝒳₀ₓ,𝒳₂ₓ = fft_integral(invx, yy, Φ′,Ψ, k, bg.ℋ(invx), q,𝕡.Σm_ν,𝕡,M)#,
        # Put it all together
        𝒳₀ = 𝒳ₛ₀ .+ real.(𝒳₀ₓ) 
        𝒳₂ = 𝒳ₛ₂ .+ real.(𝒳₂ₓ) 

        all_splines₀[i_q+2] = linear_interpolation(invx,𝒳₀)
        all_splines₂[i_q+2] = linear_interpolation(invx,𝒳₂)

    end


    return invx, all_splines₀, all_splines₂
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
                saveat=ie.bg.x_grid, dense=false, #FIXME
                )
    return sol
end

function boltsolve_flex(ie::IEallν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
    prob = ODEProblem{true}(Bolt.ie!, u₀, (x_ini , x_fin), ie)
    sol = solve(prob, ode_alg, reltol=reltol,
                saveat=ie.bg.x_grid, dense=false, #FIXME
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
function boltsolve_conformal_flex(confie::ConformalIEallν{T},#FIXME we don't need this? {Hierarchy{T},AbstractInterpolation{T}},
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
    M, reltol,x_ini, x_fin,u0,m,q,i_q) where T
    # 𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
    ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒩₀_km1, 𝒩₂_km1,
                    ℓᵧ, ℓ_mν, n_q)
    perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie(ie_k_late,perturb_k_late,M,m,q,i_q,
                        u0,perturb_k_late.t) 
    return xx,𝒩₀_k,𝒩₂_k,perturb_k_late
end

function iterate_fft_c(𝒩₀_km1,𝒩₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, ℓ_mν, n_q,
    M, reltol,η_ini, η_fin,u0,m,q,i_q) where T
    # 𝒩₀_k,𝒩₂_k = zero(𝒩₀_km1),zero(𝒩₂_km1) #need this line ow is never updated
    ie_k_late = IEν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒩₀_km1, 𝒩₂_km1,
                    ℓᵧ, ℓ_mν, n_q)
    ie_k_conf_late_c = ConformalIEν(ie_k_late,η2x);
    perturb_k_late_c = boltsolve_conformal_flex(ie_k_conf_late_c, η_ini, η_fin, u0; reltol=reltol)
    xx,𝒩₀_k,𝒩₂_k = fft_ie_c(ie_k_conf_late_c.ie,perturb_k_late_c,M,m,q,i_q,
                        u0,η2x(perturb_k_late_c.t/Mpcfac)) 
    return xx,𝒩₀_k,𝒩₂_k,perturb_k_late_c
end

function iterate_fft_allν(𝒳₀_km1,𝒳₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, n_q,
    M, reltol,x_ini, x_fin,u0) where T
    ie_k_late = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒳₀_km1,𝒳₂_km1,
                    ℓᵧ, n_q)
    perturb_k_late = boltsolve_flex(ie_k_late, x_ini, x_fin, u0; reltol=reltol)
    # ~~Loop here? 1?~~
    xx,𝒳₀_k,𝒳₂_k = fft_ie(ie_k_late,perturb_k_late,M,
                        u0,perturb_k_late.t) 
    return xx,𝒳₀_k,𝒳₂_k,perturb_k_late
end

function iterate_fft_allν_c(𝒳₀_km1,𝒳₂_km1, 𝕡::CosmoParams{T}, bg, ih, k, ℓᵧ, n_q,
    M, reltol,η_ini, η_fin,u0) where T
    ie_k_late = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
                    𝒳₀_km1,𝒳₂_km1,
                    ℓᵧ, n_q)
    ie_k_conf_late_c = ConformalIEallν(ie_k_late,η2x);
    perturb_k_late_c = boltsolve_conformal_flex(ie_k_conf_late_c, η_ini, η_fin, u0; reltol=reltol)
    xx,𝒳₀_k,𝒳₂_k = fft_ie_c(ie_k_conf_late_c.ie,perturb_k_late_c,M,
                        u0,η2x(perturb_k_late_c.t/Mpcfac)) 
    return xx,𝒳₀_k,𝒳₂_k,perturb_k_late_c
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
                                   M,reltol,x_ini,x_fin,u0,m,q,i_q)
    end
    return xx_k, 𝒩₀_k,𝒩₂_k,perturb_k
end
#ctime version
function itersolve_fft(Nₖ::Int,ie_0_c::ConformalIEν{T},M::Int,η_ini, η_fin,u0;reltol=1e-6) where T
    𝒩₀_0,𝒩₂_0 = ie_0_c.ie.s𝒩₀,ie_0_c.ie.s𝒩₂
    𝒩₀_k,𝒩₂_k = 𝒩₀_0,𝒩₂_0
    perturb_k = nothing
    ηη_k = nothing
    for k in 1:Nₖ
        ηη_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft_c(𝒩₀_k,𝒩₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.ℓ_mν,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0,m,q,i_q)
    end
    return ηη_k,𝒩₀_k,𝒩₂_k,perturb_k
end

function itersolve_fft(Nₖ::Int,ie_0::IEallν{T},M::Int,x_ini,x_fin,u0;reltol=1e-6) where T
    𝒳₀_0,𝒳₂_0 = ie_0.s𝒳₀,ie_0.s𝒳₂
    𝒳₀_k,𝒳₂_k = 𝒳₀_0,𝒳₂_0
    perturb_k = nothing
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate_fft_allν(𝒳₀_k,𝒳₂_k,ie_0.par,ie_0.bg,ie_0.ih,
                                   ie_0.k,ie_0.ℓ_γ,ie_0.nq,
                                   M,reltol,x_ini,x_fin,u0)
    end
    return xx_k,𝒳₀_k,𝒳₂_k,perturb_k
end
#ctime version
function itersolve_fft(Nₖ::Int,ie_0_c::ConformalIEallν{T},M::Int,η_ini, η_fin,u0;reltol=1e-6) where T
    𝒳₀_0,𝒳₂_0 = ie_0_c.ie.s𝒳₀,ie_0_c.ie.s𝒳₂
    𝒳₀_k,𝒳₂_k = 𝒳₀_0,𝒳₂_0
    perturb_k = nothing
    ηη_k = nothing
    for k in 1:Nₖ
        ηη_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate_fft_allν_c(𝒳₀_k,𝒳₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0)
    end
    return ηη_k,𝒳₀_k,𝒳₂_k,perturb_k
end


# Helper functon for switch
function get_switch_u0(η,hierarchy_conf) #Input is η of the switch
    # This function assumes truncated hierarchies for all neutrinos (but not yet photons)
    hierarchy,bg = hierarchy_conf.hierarchy,hierarchy_conf.hierarchy.bg
    Mpcfac = bg.H₀*299792.458/100.
    switch_idx = argmin(abs.(bg.η*Mpcfac .-η)) #for now we use the bg to find the switch
    #solve the split ode
    ℓᵧ,ℓ_ν,n_q = hierarchy.ℓᵧ,hierarchy.ℓ_ν, hierarchy.nq
    pertlen=2(ℓᵧ+1) + (ℓ_ν+1) + (ℓ_mν+1)*n_q + 5
    # \/ we want to report this timing to get a full picture of total time (early+late)
    sol_early_c = h_boltsolve_conformal_flex(hierarchy_conf, bg.η[1], bg.η[switch_idx],  initial_conditions(bg.x_grid[1], hierarchy));
    
    # Get the new initial conditions
    u0_ie_c = zeros(2(ℓᵧ+1) + (2+1) + (2+1)*n_q + 5);
    # The first split will be the same
    for i in  1:2(ℓᵧ+1)+(2+1) #up to massless ν quadrupole (w/photon hierarchy #FIXME)
        u0_ie_c[i] = sol_early_c.u[end][i]
    end

    # for i in  2(ℓᵧ+1)+(ℓ_ν+1)+1:pertlen #skip the higher massless hierarchy multipoles
    for i in  2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+1 #skip the higher massless hierarchy multipoles
        down_shift = i-(ℓ_ν-2) #shift down by the number of multipoles we skip
        u0_ie_c[down_shift] = sol_early_c.u[end][i]
    end

    #TODO Since this is contiguous, can combine it with the above loop
    # Do the same for massive neutrinos, which are arranged as [q1,q2,...,qnq]_ℓ=0, [q1,q2,...,qnq]_ℓ=1, ..., []_ℓ=ℓ_mν
    for i_ℓ in 1:3 #we fill all the multipoles as we usually would up to the quadrupole
        for i in 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*(i_ℓ-1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+i_ℓ*n_q 
            down_shift = i-(ℓ_ν-2)
            u0_ie_c[down_shift] = sol_early_c.u[end][i]
        end
    end

    for i in  pertlen-5:pertlen #skip the higher massless hierarchy multipoles
        down_shift = i-(ℓ_ν-2)-n_q*(ℓ_mν-2) #shift down by the number of multipoles we skip
        u0_ie_c[down_shift] = sol_early_c.u[end][i]
    end

    return u0_ie_c
end


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
n_q=15
bg = Background(𝕡; x_grid=ret[1,1]:round(dx,digits=3):ret[end,1], nq=n_q);
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b); #FIXME γΩ
ih = IonizationHistory(𝕣, 𝕡, bg);
Mpcfac = bg.H₀*299792.458/100.
k = Mpcfac*kclass #get k in our units


#unrelated code determining if we need "T" in all the places it is for interpolation
using ForwardDiff
ForwardDiff.Dual(0.67)
𝕡 = CosmoParams{ForwardDiff.Dual}(h=ForwardDiff.Dual(0.67));
bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=n_q);
bg

function fbg(Ω_b::DT) where DT
    𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
    println(typeof(bg))
    return bg.η(-5)
 end

 using FastGaussQuadrature
 qp,qw = gausslegendre( 15 )

 typeof(qp) 
 typeof(qp) <: Vector{T where T >: Float64}
Real >: Float64
Vector{Real} >: Vector{Float64}

typeof(qp) <: Array{Real,1}
typeof(qp) <: Array{Float64,1}

 fbg(0.046)
 Δ = 1e-3
 (fbg(0.046+ Δ) - fbg(0.046 - Δ)) / 2Δ
 ForwardDiff.derivative(fbg, 0.046)
fbg(ForwardDiff.Dual(0.046,1))

typeof(Vector{Float64}) <: typeof(Vector{Real})



# Hierarchy for comparison purposes - now replace with conformal hierarchy...
ℓᵧ=50
ℓ_mν=20
ℓ_ν=50#3#3#ℓ_ν10#ℓᵧ
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
η2x = linear_interpolation(bg.η.(bg.x_grid),bg.x_grid)
hierarchy_conf = ConformalHierarchy(hierarchy,η2x);
results_conf = boltsolve_conformal(hierarchy_conf;reltol=reltol);

#truncated conformal hierarchy
# Input to the ie integrator struct (akin to hierarchy)
𝒩₀_0,𝒩₂_0 =  results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:] #hierarchy answer
spl0h𝒩₀,spl0h𝒩₂ = linear_interpolation(bg.x_grid,𝒩₀_0), linear_interpolation(bg.x_grid,𝒩₂_0)
ν_idx = 2(ℓᵧ+1) + 1
# ie_0 = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
#         spl0h𝒩₀,
#         spl0h𝒩₂,
#         ℓᵧ, ℓ_mν, n_q);
# perturb_0 = boltsolve(ie_0;reltol=reltol); #no rsa

# ie_0_conf = ConformalIEν(ie_0,η2x);
# results_conf_ie_0 = boltsolve_conformal(ie_0_conf;reltol=reltol);

c𝒩₀_0,c𝒩₂_0 =  results_conf[2(ℓᵧ+1)+1,:],results_conf[2(ℓᵧ+1)+3,:] #hierarchy answer
c_spl0h𝒩₀,c_spl0h𝒩₂ = linear_interpolation(η2x(results_conf.t/Mpcfac),c𝒩₀_0), linear_interpolation(η2x(results_conf.t/Mpcfac),c𝒩₂_0)

cℳ₀q1_0,cℳ₂q1_0 =  results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+1,:],results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+1,:] #hierarchy answer
c_spl0hℳ₀q1,c_spl0hℳ₂q1 = linear_interpolation(η2x(results_conf.t/Mpcfac),cℳ₀q1_0), linear_interpolation(η2x(results_conf.t/Mpcfac),cℳ₂q1_0)
cℳ₀qend_0,cℳ₂qend_0 =  results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+15,:],results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+15,:] #hierarchy answer
c_spl0hℳ₀qend,c_spl0hℳ₂qend = linear_interpolation(η2x(results_conf.t/Mpcfac),cℳ₀qend_0), linear_interpolation(η2x(results_conf.t/Mpcfac),cℳ₂qend_0)

plot(bg.x_grid,c_spl0hℳ₀q1.(bg.x_grid))
plot!(bg.x_grid,c_spl0hℳ₀qend.(bg.x_grid))
plot(bg.x_grid,c_spl0hℳ₂q1.(bg.x_grid))
plot!(bg.x_grid,c_spl0hℳ₂qend.(bg.x_grid))


u0_ie_c = get_switch_u0(1.0,hierarchy_conf)

M=8192
# fft_ie_c(ie_0,perturb_0,M,𝕡.Σm_ν,q1,1,u0_ie_c,bg.x_grid);
# ℳ₀[0+i_q],ℳ₀[0+nq+i_q],ℳ₀[0+2nq+i_q]
ℓ_ν,ℓ_mν
massive_interps₀ = [linear_interpolation(η2x(results_conf.t/Mpcfac),results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]) for idx_q in 1:n_q];
massive_interps₂ = [linear_interpolation(η2x(results_conf.t/Mpcfac),results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2(ℓ_mν+1)+idx_q,:]) for idx_q in 1:n_q];
all_splines₀ = [c_spl0h𝒩₀,massive_interps₀...]
all_splines₂ = [c_spl0h𝒩₂,massive_interps₂...]
ie_all_0 = IEallν(BasicNewtonian(), 𝕡, bg, ih, k, #test the new struct
        all_splines₀,
        all_splines₂,
        ℓᵧ, n_q);
        # and test the evolution...
perturb_all_0 = boltsolve(ie_all_0;reltol=reltol);

BasicNewtonian <: Bolt.PerturbationIntegrator

typeof(ie_all_0)<: IEallν
typeof(ie_all_0)

ie_all_0_c = ConformalIEallν(ie_all_0,η2x);
perturb_all_0_c = boltsolve_conformal(ie_all_0_c;reltol=2e-5);
#^THIS IS UNSTABLE? FIXME!!! happens only for rtol<2e-5
perturb_all_0.u[:,1]

#FIMXE SOMETHING IS ALREADY WRONG NEAR ICS???
# Maybe an indexing error? Or actually physical?
plot(bg.x_grid,hcat(perturb_all_0.u...)[1,:])
plot!(η2x.(perturb_all_0_c.t/Mpcfac),hcat(perturb_all_0_c.u...)[1,:])
plot!(bg.x_grid,results[1,:])
plot(bg.x_grid,results[ν_idx,:])
plot!(bg.x_grid,all_splines₀[1].(bg.x_grid))
plot(bg.x_grid,results[ν_idx+ℓ_ν+1,:])
plot!(bg.x_grid,all_splines₀[2].(bg.x_grid))
plot(bg.x_grid,results[end-4-n_q*ℓ_mν-1,:])
plot!(bg.x_grid,all_splines₀[end].(bg.x_grid))


plot!(bg.x_grid,all_splines₀[2].(bg.x_grid))

function χ′z(a,q,m,tq_pts,tq_wts)
    return q / (a * Bolt.ℋ_a(a,𝕡,tq_pts,tq_wts) * √(q^2 + (a*m)^2 ) )
end
# WHAT WAS I DOING WITH THIS?


xx_k,𝒳₀_k,𝒳₂_k,perturb_k = itersolve_fft(5,ie_all_0_c,M,
                                            1.0/Mpcfac,bg.η[end],u0_ie_c;reltol=1e-5)

plot(bg.x_grid,𝒳₀_k[1,:])

#------------------------------------------------
#GENERALIZE THIS TO FOR LOOP OVER Q PTS
#save  neutrinos:
# writedlm("./test/data/Bolt_mslss_nuperts_nonu_lmax$(ℓ_ν).dat",
#           hcat(bg.x_grid,results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:]))
Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
logqmin,logqmax=log10(Tν/30),log10(Tν*30)
q1,q3,qmid,q10,q11,q12,qend = Bolt.xq2q(bg.quad_pts[1],logqmin,logqmax),Bolt.xq2q(bg.quad_pts[3],logqmin,logqmax),Bolt.xq2q(bg.quad_pts[8],logqmin,logqmax),Bolt.xq2q(bg.quad_pts[10],logqmin,logqmax),Bolt.xq2q(bg.quad_pts[11],logqmin,logqmax),Bolt.xq2q(bg.quad_pts[12],logqmin,logqmax),Bolt.xq2q(bg.quad_pts[end],logqmin,logqmax)
χt0 =  [Bolt.χν(x, q1 , 0.0 , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χt1 =  [Bolt.χν(x, q1 , 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χt3 =  [Bolt.χν(x, q3 , 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χtmid =  [Bolt.χν(x, qmid , 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χt10 =  [Bolt.χν(x, q10 , 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χt11 =  [Bolt.χν(x, q11, 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χt12 =  [Bolt.χν(x, q12 , 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
χtend =  [Bolt.χν(x, qend , 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
plot(bg.η*Mpcfac,bg.η*Mpcfac,ls=:dash,color=:black,label="η",legend=:topleft)
plot!(bg.η*Mpcfac,χt0*Mpcfac,label="χt0",xscale=:log10,yscale=:log10,ls=:dot)
plot!(bg.η*Mpcfac,χt1*Mpcfac,label="χt1")
plot!(bg.η*Mpcfac,χt3*Mpcfac,label="χt3")
plot!(bg.η*Mpcfac,χtmid*Mpcfac,label="χtmid")
plot!(bg.η*Mpcfac,χt10*Mpcfac,label="χt10")
plot!(bg.η*Mpcfac,χt12*Mpcfac,label="χt12")
plot!(bg.η*Mpcfac,χtend*Mpcfac,label="χtend")

yyxt1 = k.* (χt1 .- χt1[1])
dyt1=(yyxt1[end]-yyxt1[1])/(M-1)
yyt = yyxt1[1]:dyt1:yyxt1[end]
invxt1 = linear_interpolation(yyxt1,bg.x_grid).(yyt)
linear_interpolation(yyxt1,bg.x_grid)
plot!(yyxt1,sort(yyxt1))
yyxtend = k.* (χtend .- χtend[1])
dytend=(yyxtend[end]-yyxtend[1])/(M-1)
yytend = yyxtend[1]:dytend:yyxtend[end]
invxtend = linear_interpolation(yyxtend,bg.x_grid).(yytend)
plot!(yyxtend,sort(yyxtend))
sum(yyxt1[2:end]-yyxt1[1:end-1] .< 0.0)
sum((yyxt1[2:end]-yyxt1[1:end-1])[end-250:end-110] .< 0.0)

plot(bg.η*Mpcfac,yyxt1,xscale=:log10,legend=:topleft)
plot!(bg.η.(bg.x_grid[end-250:end-110]).*Mpcfac,yyxt1[end-250:end-110])
xlabel!("η")
ylabel!("kΔχ")
plot(bg.η*Mpcfac,χt1,xscale=:log10,legend=:topleft)
plot!(bg.η.(bg.x_grid[end-250:end-110]).*Mpcfac,χt1[end-250:end-110])


yyxtmid = k.* (χtmid .- χtmid[1])
dytmid=(yyxt1[end]-yyxtmid[1])/(M-1)
yytmid = yyxtmid[1]:dytmid:yyxtmid[end]
plot(yyxtmid,sort(yyxtmid))
sum(yyxtmid[2:end]-yyxtmid[1:end-1] .< 0.0)


yyxt3 = k.* (χt3.- χt3[1])
dyt3=(yyxt3[end]-yyxt3[1])/(M-1)
yyt3= yyxt3[1]:dyt3:yyxt3[end]
plot(yyxt3,sort(yyxt3))
sum(yyxt3[2:end]-yyxt3[1:end-1] .< 0.0)
linear_interpolation(yyxt3,bg.x_grid)


yyxt10 = k.* (χt10.- χt10[1])
dyt10=(yyxt10[end]-yyxt10[1])/(M-1)
yyt10= yyxt10[1]:dyt10:yyxt10[end]
sum(yyxt10[2:end]-yyxt10[1:end-1] .< 0.0)
linear_interpolation(yyxt10,bg.x_grid)

# 11 is the first q for which the neutrino horizon is actually monotonic
yyxt11 = k.* (χt11.- χt11[1])
dyt11=(yyxt11[end]-yyxt11[1])/(M-1)
yyt11= yyxt11[1]:dyt11:yyxt11[end]
sum(yyxt11[2:end]-yyxt11[1:end-1] .< 0.0)
linear_interpolation(yyxt11,bg.x_grid)

yyxt12 = k.* (χt12.- χt12[1])
dyt12=(yyxt12[end]-yyxt12[1])/(M-1)
yyt12= yyxt12[1]:dyt12:yyxt12[end]
sum(yyxt12[2:end]-yyxt12[1:end-1] .< 0.0)
linear_interpolation(yyxt12,bg.x_grid)


#Messing with neutrino horizon
function χν_old(x, q, m, par::AbstractCosmoParams,quad_pts,quad_wts) 
    # adding m here is a bit annoying but we need the ability to use massless neutrinos
    logamin,logamax=-13.75,log10(Bolt.x2a(x)) #0,x2a(x)
    ϵ(a,q) = √(q^2 + (a*m)^2 )
    Iχν(y) = 1.0 / (Bolt.xq2q(y,logamin,logamax) * Bolt.ℋ_a(Bolt.xq2q(y,logamin,logamax), par,quad_pts,quad_wts) * ϵ(Bolt.xq2q(y,logamin,logamax),q)
                   )/ Bolt.dxdq(Bolt.xq2q(y,logamin,logamax),logamin,logamax)
    return q*sum(Iχν.(quad_pts).*quad_wts)
end

χν_old(bg.x_grid[end], q1, 𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts)

𝕡.Σm_ν/q1,𝕡.Σm_ν/qend

#Let's schematically look at the integrand for the final x since this is around where there are issues
plot(bg.η.*Mpcfac,q1./(q1^2 .+ (Bolt.x2a.(bg.x_grid)*𝕡.Σm_ν).^2 ).^(1/2),xscale=:log10)
plot!(bg.η.*Mpcfac,qend./(qend^2 .+ (Bolt.x2a.(bg.x_grid)*𝕡.Σm_ν).^2 ).^(1/2),xscale=:log10)

#now do it the way we do it in log10a
final_idx = length(bg.x_grid) #- 175 #260 #pick some intermediate index
logamin,logamax=-13.75,bg.x_grid[final_idx]/log(10.)#0.0
plot(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .* q1./(q1^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
    #  label=L"$\chi(q_{i})$ low (260)",#yscale=:log10,
     label=L"$\chi(q_{i})$ high (175)",ls=:dash,#yscale=:log10,
     legend=:topleft,left_margin=4mm)
plot!(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*q3./(q3^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{3})$" )
plot!(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*qmid./(qmid^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{8})$" )
plot!(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*q10./(q10^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{10})$" )
plot!(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*q11./(q11^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{11})$" )
plot!(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*qend./(qend^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{f})$" )
plot!(bg.x_grid[1:final_idx]./log(10.),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),label=L"$\eta$" )
quadpts_log10a = log10.(Bolt.xq2q.(bg.quad_pts,logamin,logamax)) #1.0 ./Bolt.dxdq.(Bolt.xq2q.(bg.quad_pts,logamin,logamax),logamin,logamax)
vline!(quadpts_log10a,
        # label="quad pts low (260)",
        label="quad pts high (175)",ls=:dash
        )
vline!([bg.x_grid[end-250]./log(10.),bg.x_grid[end-110]./log(10.)],color=:red,label="q1 problem zone")
xlabel!(L"$\log_{10}(a)$")
ylabel!(L"$\eta_{f}$ ctime integrand")
xlims!(-2.0,0.0)
plot!(legend=:topright)
plot!(yscale=:log10)
ylims!(3e30,5e31)
savefig("../misc_plots/fft_debug/fft_experiments/mssv_chi_q1_integrand_zoom_log_zoom.pdf")
bg.quad_wts

plot(bg.x_grid,bg.ℋ,yscale=:log10)


q1/q3

# Think about a change of variable, because then we only do this integral once
plot!(bg.x_grid[1:final_idx]./log(10.) .+ log10.(𝕡.Σm_ν/q1) ,(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .* q1./(q1^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{i})$ high (175)",ls=:dash,yscale=:log10,
     legend=:topleft,left_margin=4mm)
plot!(bg.x_grid[1:final_idx]./log(10.) .+ log10.(𝕡.Σm_ν/q3), (1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*q3./(q3^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{3})$" )
plot!(bg.x_grid[1:final_idx]./log(10.) .+ log10.(𝕡.Σm_ν/qmid),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*qmid./(qmid^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{8})$" )
plot!(bg.x_grid[1:final_idx]./log(10.) .+ log10.(𝕡.Σm_ν/q10),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*q10./(q10^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{10})$" )
plot!(bg.x_grid[1:final_idx]./log(10.) .+ log10.(𝕡.Σm_ν/q11),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*q11./(q11^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{11})$" )
plot!(bg.x_grid[1:final_idx]./log(10.) .+ log10.(𝕡.Σm_ν/qend),(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) .*qend./(qend^2 .+ (Bolt.x2a.(bg.x_grid[1:final_idx])*𝕡.Σm_ν).^2 ).^(1/2) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),
     label=L"$\chi(q_{f})$" )
plot!(bg.x_grid[1:final_idx]./log(10.) ,(1.0 ./ (Bolt.x2a.(bg.x_grid[1:final_idx]).*bg.ℋ.(bg.x_grid[1:final_idx])) ) ./ Bolt.dxdq.(exp.(bg.x_grid[1:final_idx]),logamin,logamax),label=L"$\eta$" )

# Ok so we can do GaussHermite in log z for a gaussian with mean zero with interval somehow rescaled.
using FastGaussQuadrature
tq_pts, tq_wts =  gausslobatto( n_q )
tq_wts

logzmin_1, logzmax_1 = logamin+ log10.(𝕡.Σm_ν/q1), logamax + log10.(𝕡.Σm_ν/q1)
ϵ(a,q) = √(q^2 + (a*𝕡.Σm_ν)^2 )
Itq(y) = 1.0 / (Bolt.xq2q(y,logzmin_1,logzmax_1) * Bolt.ℋ_a(
                Bolt.xq2q(y,logzmin_1,logzmax_1), 𝕡,tq_pts, tq_wts
                    ) * ϵ(Bolt.xq2q(y,logzmin_1,logzmax_1),q1)
               )/ Bolt.dxdq(Bolt.xq2q(y,logzmin_1,logzmax_1),logzmin_1,logzmax_1)

q1*sum(Itq.(tq_pts).*tq_wts)
q1*sum(Itq.(bg.quad_pts).*bg.quad_wts)
χt1[final_idx]

log10.(Bolt.xq2q.(tq_pts,logzmin_1,logzmax_1)) #takes unit interval to logz
scatter(1:1:length(tq_pts),log10.(Bolt.xq2q.(tq_pts,logzmin_1,logzmax_1)) )
ylims!(-35,35)
hline!([Bolt.xq2q.(tq_pts,logzmin_1,logzmax_1)[8]])

function χν_new(x, q, m, par::AbstractCosmoParams,quad_pts,quad_wts) 
    # adding m here is a bit annoying but we need the ability to use massless neutrinos
    logamin,logamax=-17.75,log10(Bolt.x2a(x)) #0,x2a(x)
    logzmin_1, logzmax_1 = logamin+ log10(m/q), logamax + log10(m/q)
    ϵ(a,q) = √(q^2 + (a*m)^2 )
    Itq(y) = 1.0 / (Bolt.xq2q(y,logzmin_1,logzmax_1) * Bolt.ℋ_a(
                    Bolt.xq2q(y,logamin,logamax), par,quad_pts, quad_wts
                        ) * ϵ(Bolt.xq2q(y,logamin,logamax),q)
                )/ Bolt.dxdq(Bolt.xq2q(y,logzmin_1,logzmax_1),logzmin_1,logzmax_1)
    return q*sum(Itq.(quad_pts).*quad_wts)
end

function χ′z(a,q,m)
    return q / (a * Bolt.ℋ_a(a,𝕡,tq_pts,tq_wts) * √(q^2 + (a*m)^2 ) )
end

log10.(𝕡.Σm_ν/q1)

χthermite = [χν_new(x, q1 , 𝕡.Σm_ν , 𝕡 ,tq_pts,tq_wts) for x in bg.x_grid]

plot(bg.η*Mpcfac,χt1*Mpcfac,label="χt1",xscale=:log10,yscale=:log10,legend=:bottomright)
plot!(bg.η*Mpcfac,χthermite*Mpcfac,label="χt1",xscale=:log10,yscale=:log10)
ylims!(10,2e2)
xlims!(1e2,2e4)

zz = 10.0.^(-4.0:0.01:0.0)
log(zz[1]*q1/(𝕡.Σm_ν))
plot(zz,1.0./zz .* 1.0./sqrt.(1.0.+zz)./Bolt.ℋ_a.(zz./(𝕡.Σm_ν/q1),𝕡,tq_pts,tq_wts) .*( (𝕡.Σm_ν/q1) ./zz)  ,
xscale=:log10,yscale=:log10)# 

zz./(𝕡.Σm_ν/q1)
# Let us do the integral in z once and for all as a cumsum type thing
# Then we can just interpolate that at the requisite values of q etc.
#no you can't do this because the integral is dz 1/z 1/sqrt(1+z) 1/H(a), but the cnxn
#btwn a and z is depedent on q/m


Bolt.ℋ_a.(zz./(𝕡.Σm_ν/q1),𝕡,tq_pts,tq_wts)


#Why are we evenn bothering with quadrature??
#Why do we not just do something like we did for optical depth??
# τ_primes = [τ′(x_, Xₑ_function, par, ℋ_function) for x_ in x]
# τ_integrated = reverse(cumul_integrate(rx, reverse(τ_primes)))
# τ̂ = interpolate((x,),τ_integrated,Gridded(Linear()))
using NumericalIntegration
χ_primes = [χ′z(a,q1,𝕡.Σm_ν) for a in aa]
aa = 10.0.^(-13.0:0.01:0.0)
χ_integrated = cumul_integrate(aa, χ_primes)
@btime χ_integrated_x = cumul_integrate(exp.(bg.x_grid),  [χ′z(exp(x),q1,𝕡.Σm_ν) for x in bg.x_grid])
#  4.227 ms (26029 allocations: 940.72 KiB)
@btime χνs = [Bolt.χν(x, q1,𝕡.Σm_ν , 𝕡 ,bg.quad_pts,bg.quad_wts) for x in bg.x_grid]
# 63.653 ms (46029 allocations: 8.05 MiB) #how is it possible this takes longer??
# well duh, 2000*12 bs 2000*1

# I am not sure there is a way to do what I want to do without specifying q,a
# We eventually have to loop over all the qs anyways (I THINKN???)
# so this is probably a misguided attempt at savings
plot(log.(aa)[2:end],χ_integrated[2:end]*Mpcfac,yscale=:log10,legend=:bottomright)
plot!(bg.x_grid[2:end],χ_integrated_x[2:end]*Mpcfac,yscale=:log10,legend=:bottomright)

plot!(bg.x_grid,χt1*Mpcfac,label="χt0",yscale=:log10,ls=:dash)
xlims!(-5.0,0.0)
ylims!(10,3e2)
sum( ((χ_integrated[2:end]*Mpcfac)[2:end]-(χ_integrated[2:end]*Mpcfac)[1:end-1]) .<0.0 )
sum( ((χ_integrated_x[2:end]*Mpcfac)[2:end]-(χ_integrated_x[2:end]*Mpcfac)[1:end-1]) .<0.0 )
#looks good
