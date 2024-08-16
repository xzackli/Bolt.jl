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
using NumericalIntegration


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
    # 𝒩[0] = ie.s𝒳₀[1](x)
    # 𝒩[2] = ie.s𝒳₂[1](x)
    𝒩₀ = ie.s𝒳₀[1](x)
    𝒩₂ = ie.s𝒳₂[1](x)
    ℳ₀ = zeros(T,nq)
    ℳ₂ = zeros(T,nq)
    for idx_q in 1:nq#0:(nq-1)
        # ℳ[0*nq+idx_q] = ie.s𝒳₀[idx_q+2](x)
        # ℳ[2*nq+idx_q] = ie.s𝒳₂[idx_q+2](x)
        ℳ₀[idx_q] = ie.s𝒳₀[idx_q+1](x)
        ℳ₂[idx_q] = ie.s𝒳₂[idx_q+1](x)
    end
    # ρℳ, σℳ  =  @views ρ_σ(ℳ[0:nq-1], ℳ[2*nq:3*nq-1], bg, a, par) #monopole (energy density, 00 part),quadrupole (shear stress, ij part)
    ρℳ, σℳ  =  @views ρ_σ(ℳ₀, ℳ₂, bg, a, par)
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
    #empty splines
    all_splines₀ = copy(ie.s𝒳₀)
    all_splines₂ = copy(ie.s𝒳₂)
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

    𝒩₀ = ie.s𝒳₀[1](x_grid[1])
    _,_,𝒩₁,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒩₂ = ie.s𝒳₂[1](x_grid[1])
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩₁,𝒩₂)) #massless
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
        ℳ₀ = ie.s𝒳₀[2+i_q](x_grid[1])
        _,_,_, ℳ₁,_,_,_,_,_ =  unpack(u₀,ie)   
        ℳ₂ = ie.s𝒳₂[2+i_q](x_grid[1])

        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀,ℳ₁[1+i_q],ℳ₂)) #massive
        
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

function fft_ie_c(ie::IEallν,perturb,M,u₀,x_grid) #FIXME add type decorators
    𝕡,bg,k,nq = ie.par,ie.bg,ie.k,ie.nq
    Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *Bolt.ρ_crit(𝕡) *𝕡.Ω_r)^(1/4)
    logqmin,logqmax=log10(Tν/30),log10(Tν*30)
    all_splines₀ = copy(ie.s𝒳₀)
    all_splines₂ = copy(ie.s𝒳₂)
    #explicitly do massless case
    χνs = cumul_integrate(exp.(x_grid),  [χ′z(exp(x),1.0,0.0,bg.quad_pts,bg.quad_wts) for x in x_grid])
    yyx = k.* (χνs .- χνs[1])
    dy=(yyx[end]-yyx[1])/(M-1)
    yy = yyx[1]:dy:yyx[end]
    invx = linear_interpolation(yyx,x_grid).(yy) #get xpoints at equispaced "neutrino ctime" FIXME use spline?
    # Get metric sources #FIXME this should probably happen outside of this function
    Φ′,Ψ = zeros(M),zeros(M)
    for j in 1:M
        Φ′[j],Ψ[j] = get_Φ′_Ψ(perturb( bg.η(invx[j]) .*Mpcfac ),ie,invx[j])
    end
    
    𝒩₀ = ie.s𝒳₀[1](x_grid[1])
    _,_,𝒩₁,_,_,_,_,_,_ =  unpack(u₀,ie)   
    𝒩₂ = ie.s𝒳₂[1](x_grid[1])
    𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,𝒩₀,𝒩₁,𝒩₂)) #massless
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
        ℳ₀ = ie.s𝒳₀[2+i_q](x_grid[1])
        _,_,_, ℳ₁,_,_,_,_,_ =  unpack(u₀,ie)   
        ℳ₂ = ie.s𝒳₂[2+i_q](x_grid[1])
        𝒳ₛ₀, 𝒳ₛ₂ = unzip(Wsum.(yy,ℳ₀,ℳ₁[1+i_q],ℳ₂)) #massive
        
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
                # saveat=ie.bg.x_grid, 
                dense=false, #FIXME
                )
    return sol
end

function boltsolve_flex(ie::IEallν{T}, x_ini,x_fin, u₀, ode_alg=KenCarp4(); reltol=1e-6) where T 
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
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒩₀_k,𝒩₂_k,perturb_k = iterate_fft_c(𝒩₀_k,𝒩₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.ℓ_mν,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0,m,q,i_q)
    end
    return xx_k,𝒩₀_k,𝒩₂_k,perturb_k
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
    xx_k = nothing
    for k in 1:Nₖ
        xx_k,𝒳₀_k,𝒳₂_k,perturb_k = iterate_fft_allν_c(𝒳₀_k,𝒳₂_k,ie_0_c.ie.par,ie_0_c.ie.bg,ie_0_c.ie.ih,
                                               ie_0_c.ie.k,ie_0_c.ie.ℓ_γ,ie_0_c.ie.nq,M,reltol,
                                               η_ini, η_fin,u0)
    end
    return xx_k,𝒳₀_k,𝒳₂_k,perturb_k
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
    u0_ie_c = zeros(2(ℓᵧ+1) + (0+1) + (0+1)*n_q + 5);
    # The first split will be the same
    for i in  1:2(ℓᵧ+1) #up to massless ν quadrupole (w/photon hierarchy #FIXME)
        u0_ie_c[i] = sol_early_c.u[end][i]
    end
    #set the massless neutrino dipole
    u0_ie_c[2(ℓᵧ+1)+1] = sol_early_c.u[end][2(ℓᵧ+1)+2]

    #This does nothing??
    # for i in  2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+1 #skip the higher massless hierarchy multipoles
    #     down_shift = i-(ℓ_ν-2) #shift the hierarchy index down by the number of multipoles we skip in the ie
    #     u0_ie_c[down_shift] = sol_early_c.u[end][i]
    # end

    #massive neutrinos, now we just do the dipole again

    #TODO Since this is contiguous, can combine it with the above loop
    # Do the same for massive neutrinos, which are arranged as [q1,q2,...,qnq]_ℓ=0, [q1,q2,...,qnq]_ℓ=1, ..., []_ℓ=ℓ_mν
    # for i_ℓ in 1:3 #we fill all the multipoles as we usually would up to the quadrupole
    #     for i in 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*(i_ℓ-1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+i_ℓ*n_q 
    #         down_shift = i-(ℓ_ν-2)
    #         u0_ie_c[down_shift] = sol_early_c.u[end][i]
    #     end
    # end

    # start at the dipole first q idx, go up to the last dipole q idx (in the hierarchy)   
    for i in 1:n_q 
        # down_shift = i-(ℓ_ν+1 - 1)-(ℓ_mν+1 - 1)
        u0_ie_c[2(ℓᵧ+1)+1+i] = sol_early_c.u[end][2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1+i]
    end

    for i in  1:5 #skip the higher massless hierarchy multipoles
        # down_shift = i-(ℓ_ν+1-1)-n_q*(ℓ_mν+1-1) #shift down by the number of multipoles we skip
        u0_ie_c[2(ℓᵧ+1)+1+n_q+i] = sol_early_c.u[end][pertlen-5+i]
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
η2x = linear_interpolation(bg.η,bg.x_grid);

hierarchy_conf = ConformalHierarchy(hierarchy,η2x);
results_conf = boltsolve_conformal(hierarchy_conf;reltol=reltol);

plot(bg.x_grid,results[2(ℓᵧ+1)+(ℓ_ν+1)+1,:])
results[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]

# writedlm("./test/data/Bolt_mslss_mssv_nuperts_mnu0p15_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat",
#           hcat(bg.x_grid,results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:],
#           [results[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:] for idx_q in 1:n_q]...,
#           [results[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+idx_q,:] for idx_q in 1:n_q]...))

# fullplanck_ansatz = readdlm("./test/data/Bolt_mslss_mssv_nuperts_nonu_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat")
fullplanck_ansatz = readdlm("./test/data/Bolt_mslss_mssv_nuperts_mnu0p15_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat")

plot(fullplanck_ansatz[:,1],fullplanck_ansatz[:,2])
plot!(fullplanck_ansatz[:,1],x_all_splines₀[1].(fullplanck_ansatz[:,1]))

plot(fullplanck_ansatz[:,1],fullplanck_ansatz[:,1+2+1])
plot!(fullplanck_ansatz[:,1],x_all_splines₀[2].(fullplanck_ansatz[:,1]))

plot(fullplanck_ansatz[:,1],fullplanck_ansatz[:,1+2+15])
plot!(fullplanck_ansatz[:,1],x_all_splines₀[16].(fullplanck_ansatz[:,1]))

plot(fullplanck_ansatz[:,1],fullplanck_ansatz[:,3])
plot!(fullplanck_ansatz[:,1],x_all_splines₂[1].(fullplanck_ansatz[:,1]))

plot(fullplanck_ansatz[:,1],fullplanck_ansatz[:,1+2+15+1])
plot!(fullplanck_ansatz[:,1],x_all_splines₂[2].(fullplanck_ansatz[:,1]))

plot(fullplanck_ansatz[:,1],fullplanck_ansatz[:,1+2+15+15])
plot!(fullplanck_ansatz[:,1],x_all_splines₂[16].(fullplanck_ansatz[:,1]))



#sometimes this happens at the end, sometimes at the beginnning....
plot(η2x.(results_conf.t/Mpcfac),results_conf(results_conf.t)[2(ℓᵧ+1)+1,:])
plot!(bg.x_grid,results[2(ℓᵧ+1)+1,:])

#truncated conformal hierarchy
𝒩₀_0,𝒩₂_0 =  results[2(ℓᵧ+1)+1,:],results[2(ℓᵧ+1)+3,:] #hierarchy answer
spl0h𝒩₀,spl0h𝒩₂ = linear_interpolation(collect(bg.x_grid),𝒩₀_0), linear_interpolation(collect(bg.x_grid),𝒩₂_0)
ν_idx = 2(ℓᵧ+1) + 1
c𝒩₀_0,c𝒩₂_0 =  results_conf[2(ℓᵧ+1)+1,:],results_conf[2(ℓᵧ+1)+3,:] #hierarchy answer
c_spl0h𝒩₀,c_spl0h𝒩₂ = linear_interpolation(η2x(results_conf.t/Mpcfac),c𝒩₀_0), linear_interpolation(η2x(results_conf.t/Mpcfac),c𝒩₂_0)

massive_interps₀ = [linear_interpolation(η2x(results_conf.t/Mpcfac),results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]) for idx_q in 1:n_q];
massive_interps₂ = [linear_interpolation(η2x(results_conf.t/Mpcfac),results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+idx_q,:]) for idx_q in 1:n_q];
all_splines₀ = [c_spl0h𝒩₀,massive_interps₀...]
all_splines₂ = [c_spl0h𝒩₂,massive_interps₂...]

x_massive_interps₀ = [linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+(ℓ_ν+1)+idx_q,:]) for idx_q in 1:n_q];
x_massive_interps₂ = [linear_interpolation(collect(bg.x_grid),results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+idx_q,:]) for idx_q in 1:n_q];
x_all_splines₀ = [spl0h𝒩₀,x_massive_interps₀...]
x_all_splines₂ = [spl0h𝒩₂,x_massive_interps₂...]

plot(bg.x_grid,c_spl0h𝒩₀.(bg.x_grid))
plot!(bg.x_grid,spl0h𝒩₀.(bg.x_grid))

ie_all_0 = IEallν(BasicNewtonian(), 𝕡, bg, ih, k, #test the new struct
        all_splines₀,
        all_splines₂,
        ℓᵧ, n_q);
        # and test the evolution...
perturb_all_0 = boltsolve(ie_all_0;reltol=reltol);
x_ie_all_0 = IEallν(BasicNewtonian(), 𝕡, bg, ih, k, #test the new struct
        x_all_splines₀,
        x_all_splines₂,
        ℓᵧ, n_q);
x_perturb_all_0 = boltsolve(x_ie_all_0;reltol=reltol);

ie_all_0_c = ConformalIEallν(ie_all_0,η2x);
perturb_all_0_c = boltsolve_conformal(ie_all_0_c;reltol=reltol);


x_ie_all_0_c = ConformalIEallν(x_ie_all_0,η2x);
x_perturb_all_0_c = boltsolve_conformal(x_ie_all_0_c;reltol=reltol);
#^This line still gives dtmin error on first solve...

function χ′z(a,q,m,tq_pts,tq_wts)
    return q / (a * Bolt.ℋ_a(a,𝕡,tq_pts,tq_wts) * √(q^2 + (a*m)^2 ) )
end
M=2048*4


#test
u0_iter = Bolt.initial_conditions(-20.0, x_ie_all_0);
xx_k,𝒳₀_k,𝒳₂_k,perturb_k = itersolve_fft(1,x_ie_all_0,M, 
                                            bg.x_grid[1],bg.x_grid[end],
                                            u0_iter
                                            );
 

plot(xx_k,𝒳₀_k[1].(xx_k))
plot!(xx_k,x_ie_all_0.s𝒳₀[1].(xx_k))

u0_iter_c = Bolt.initial_conditions(η2x(1.0/Mpcfac), ie_all_0_c.ie);


#something is wrong with this...it should not be going to x=-20
xx_k_c,𝒳₀_k_c,𝒳₂_k_c,perturb_k_c = itersolve_fft(2,ie_all_0_c,M, 
                                            1.0/Mpcfac,bg.η[end],
                                            u0_iter_c
                                            );
plot(xx_k_c,𝒳₀_k_c[1].(xx_k_c))
plot!(xx_k_c,ie_all_0_c.ie.s𝒳₀[1].(xx_k_c))
               

#^So neither of these things work, which is maybe not suprising since the
#ICs that have to be set at x=-20 are being set at x=-12.3 so we should maybe
#expect problems...
#Not sure why this gets so much worse at second iter...

#---------------------------------#
#---------------------------------#
# Begin Experiments
#---------------------------------#
#---------------------------------#


u0_ie_c = get_switch_u0(1.0,hierarchy_conf);
# length(u0_ie_c)-5-15
# ν_idx
# Set up ansatzs and struct
# don't forget that u0 no longer has monopole and quadrupole for truncated heirarchy, so use hierarchy results...
all_const_ansatz₀ =[linear_interpolation(η2x.(results_conf.t/Mpcfac),results_conf[ν_idx,1]*ones(length(results_conf.t))), 
                    [linear_interpolation(η2x.(results_conf.t/Mpcfac),results_conf[ν_idx+(ℓ_ν+1)+(idx_q-1),1]*ones(length(results_conf.t)))
                 for idx_q in 1:n_q]...];
# all_zero_ansatz₂ = [linear_interpolation(η2x.(results_conf.t/Mpcfac),zeros(length(results_conf.t)))
#                 for idx_q in 1:n_q+1];
all_const_ansatz₂ =[linear_interpolation(η2x.(results_conf.t/Mpcfac),results_conf[ν_idx+2,1]*ones(length(results_conf.t))), 
                    [linear_interpolation(η2x.(results_conf.t/Mpcfac),results_conf[ν_idx+(ℓ_ν+1)+2n_q+(idx_q-1),1]*ones(length(results_conf.t)))
                 for idx_q in 1:n_q]...];

x_all_const_ansatz₀ =[linear_interpolation(collect(bg.x_grid),results[ν_idx,1]*ones(length(collect(bg.x_grid)))), 
                 [linear_interpolation(collect(bg.x_grid),results[ν_idx+(ℓ_ν+1)+(idx_q-1),1]*ones(length(collect(bg.x_grid))))
              for idx_q in 1:n_q]...];

x_all_const_ansatz₂ =[linear_interpolation(collect(bg.x_grid),              results[ν_idx+2,1]*ones(length(collect(bg.x_grid)))), 
                 [linear_interpolation(collect(bg.x_grid),results[ν_idx+(ℓ_ν+1)+2n_q+(idx_q-1),1]*ones(length(collect(bg.x_grid))))
              for idx_q in 1:n_q]...];

only_massless_ansatz₀ =[
                    linear_interpolation(η2x.(results_conf.t/Mpcfac),results_conf[ν_idx,1]*ones(length(results_conf.t))), 
                    massive_interps₀...];
only_massless_ansatz₂ =[
                    linear_interpolation(η2x.(results_conf.t/Mpcfac),results_conf[ν_idx+2,1]*ones(length(results_conf.t))), 
                    massive_interps₂...];

x_only_massless_ansatz₀ =[
                        linear_interpolation(collect(bg.x_grid),results_conf[ν_idx,1]*ones(length(collect(bg.x_grid)))), 
                        x_massive_interps₀...];
x_only_massless_ansatz₂ =[
                        linear_interpolation(collect(bg.x_grid),results[ν_idx+2,1]*ones(length(collect(bg.x_grid)))), 
                        x_massive_interps₂...];
    
ie_0_late_c = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
                    # only_massless_ansatz₀,
                    # only_massless_ansatz₂,
                    # all_const_ansatz₀,
                    # all_zero_ansatz₂,
                    # all_const_ansatz₂,
                    all_splines₀,
                    all_splines₂,
                    ℓᵧ, n_q);
ie_0_conf_late_c = ConformalIEallν(ie_0_late_c,η2x);



#planck ansatzs 
planck_nonu_ansatz_data = readdlm("./test/data/Bolt_mslss_mssv_nuperts_nonu_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat");
planck_heavynu_ansatz_data = readdlm("./test/data/Bolt_mslss_mssv_nuperts_mnu0p15_msslsslmax$(ℓ_ν)_mssvlmax$(ℓ_mν).dat"); #0.15eV

planck_nonu_ansatz₀ = [linear_interpolation(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz_data[:,2]),
                        [linear_interpolation(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz_data[:,1+2+idx_q]) for idx_q in 1:n_q]...];

planck_nonu_ansatz₂ = [linear_interpolation(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz_data[:,3]),
                        [linear_interpolation(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz_data[:,1+2+n_q+idx_q]) for idx_q in 1:n_q]...];

planck_heavynu_ansatz₀ = [linear_interpolation(planck_nonu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,2]),
                        [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,1+2+idx_q]) for idx_q in 1:n_q]...];

planck_heavynu_ansatz₂ =  [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,3]),
                        [linear_interpolation(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz_data[:,1+2+n_q+idx_q]) for idx_q in 1:n_q]...];         


plot(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz₀[1].(planck_nonu_ansatz_data[:,1]))
plot!(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz₀[1].(planck_heavynu_ansatz_data[:,1]))

plot(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz₀[3].(planck_nonu_ansatz_data[:,1]))
plot!(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz₀[3].(planck_heavynu_ansatz_data[:,1]))


plot(planck_nonu_ansatz_data[:,1],planck_nonu_ansatz₀[16].(planck_nonu_ansatz_data[:,1]))
plot!(planck_heavynu_ansatz_data[:,1],planck_heavynu_ansatz₀[16].(planck_heavynu_ansatz_data[:,1]))


x_ie_0_late_c = IEallν(BasicNewtonian(), 𝕡, bg, ih, k,
                    # planck_nonu_ansatz₀,
                    # planck_nonu_ansatz₂,
                    planck_heavynu_ansatz₀,
                    planck_heavynu_ansatz₂,
                    # x_all_const_ansatz₀,
                    # x_all_const_ansatz₂,
                    # x_only_massless_ansatz₀,
                    # x_only_massless_ansatz₂,
                    # x_all_splines₀,
                    # x_all_splines₂,
                    ℓᵧ, n_q);
x_ie_0_conf_late_c = ConformalIEallν(x_ie_0_late_c,η2x);


# Experiments


# plot(xx_kt,all_const_ansatz₀[1].(xx_kt))
plot(xx_k,spl0h𝒩₀.(xx_k))
plot!(xx_k,c_spl0h𝒩₀.(xx_k))
plot!(xx_k,only_massless_ansatz₀[1].(xx_k))
plot!(xx_k,all_splines₀[1].(xx_k))

plot(xx_k,spl0h𝒩₂.(xx_k))
plot!(xx_k,c_spl0h𝒩₂.(xx_k))
plot!(xx_k,only_massless_ansatz₂[1].(xx_k))
plot!(xx_k,all_splines₂[1].(xx_k))

plot(xx_k,massive_interps₀[end].(xx_k))
plot!(xx_k,x_massive_interps₀[end].(xx_k))
plot!(xx_k,only_massless_ansatz₀[end].(xx_k))

plot(xx_k,spl0h𝒩₂.(xx_k))
plot!(xx_k,c_spl0h𝒩₂.(xx_k))
plot!(xx_k,only_massless_ansatz₂[1].(xx_k))


plot(xx_kt,all_const_ansatz₂[1].(xx_kt))
plot!(xx_kt,spl0h𝒩₂.(xx_kt))
plot(xx_kt,abs.(all_const_ansatz₀[2].(xx_kt)),yscale=:log10)
plot!(xx_kt,abs.(all_splines₀[2].(xx_kt)))
ylims!(1e-3,5e-2)
plot(xx_kt,abs.(all_const_ansatz₂[2].(xx_kt)),yscale=:log10)
plot!(xx_kt,abs.(all_splines₂[2].(xx_kt)))
plot!(xx_kt,abs.(x_all_splines₂[2].(xx_kt)))
plot!(xx_kt,abs.(results_conf(bg.η(xx_kt)*Mpcfac)[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+1,:]))
plot!(bg.x_grid,abs.(results_conf(bg.η(bg.x_grid)*Mpcfac)[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+1,:]))

# results_conf[2(ℓᵧ+1)+(ℓ_ν+1)+2n_q+1,2]


plot(xx_kt,abs.(all_const_ansatz₂[end].(xx_kt)),yscale=:log10)
plot!(xx_kt,abs.(all_splines₂[end-1].(xx_kt)))


plot(xx_kt,abs.(all_const_ansatz₂[end].(xx_kt)),yscale=:log10)
plot!(xx_kt,abs.(all_splines₂[end].(xx_kt)))

all_const_ansatz₀[1].(xx_kt)
results[end,1]
u0_ie_c[end]


results[103+1,1]
u0_ie_c[103]

results[103+50+21,1]
u0_ie_c[104]

#Ah man u0_ie_c is messed up...


# And now ctime
M=2048*4
# Nᵢ=1
reltol=1e-12
#changing k, switch, hierarchy truncation, and ansatz will need to have a re-doing of u0_ie, ie_0 struct

# First experiment
η_switchη_switch = [5.0] #[0.5,1.0,10.0,100.0] 
η_switch = 5.0
MM = [8192]#[2^i for i in 12:14]
NᵢNᵢ = [1,5]#[2i-1 for i in 1:5] #max iters

#run this for plotting consistency
xx_kt,𝒳₀_kt,𝒳₂_kt,perturb_kt= itersolve_fft(1,ie_0_conf_late_c,MM[end],
    η_switchη_switch[1]/Mpcfac,bg.η[end],get_switch_u0(η_switchη_switch[1],hierarchy_conf);reltol=reltol);

plot_idx=16;
# for η_switch in η_switchη_switch
    # Set the initial conditions at a particular switch value

    u0_ie_c = get_switch_u0(η_switch,hierarchy_conf)
    # Hierarchy
    # p1 = plot(xx_kt,all_splines₀[plot_idx].(xx_kt),label="H",color=:black,lw=2)
    # p2 = plot(xx_kt,all_splines₂[plot_idx].(xx_kt),label=false,color=:black,lw=2)
    p1 = plot(xx_kt,abs.(all_splines₀[plot_idx].(xx_kt)),label="H",color=:black,lw=2,yscale=:log10)
    p2 = plot(xx_kt,abs.(all_splines₂[plot_idx].(xx_kt)),label=false,color=:black,lw=2,yscale=:log10)

    # Initial guess (for now jsut leave the monopole)
    plot!(p1,xx_kt,1e-16 .+ abs.(x_ie_0_late_c.s𝒳₀[plot_idx].(xx_kt)),label="I = 0, mono init ansatz",legendfont=font(4),ls=:dash,yscale=:log10)
    plot!(p2,xx_kt,1e-16 .+abs.(x_ie_0_late_c.s𝒳₂[plot_idx].(xx_kt)),label=false,c=p1.series_list[2][:linecolor],legendfont=font(4),ls=:dash,yscale=:log10)
    # plot!(p1,xx_kt,x_ie_0_late_c.s𝒳₀[plot_idx].(xx_kt),label="I = 0, mono init ansatz",legendfont=font(4),ls=:dash)
    # plot!(p2,xx_kt,x_ie_0_late_c.s𝒳₂[plot_idx].(xx_kt),label=false,c=p1.series_list[2][:linecolor],legendfont=font(4),ls=:dash)
 

    
    for M in MM
        for Nᵢ in NᵢNᵢ
            xx_k,𝒳₀_k,𝒳₂_k,perturb_k = itersolve_fft(Nᵢ,
                                            # x_ie_0_conf_late_c,
                                            x_ie_0_late_c,
                                            M,
                                            # η_switch/Mpcfac,bg.η[end],
                                            η2x(η_switch/Mpcfac),bg.x_grid[end],
                                            u0_ie_c;reltol=reltol);

            println("(M = $M, Nᵢ = $Nᵢ), 
                    error against full hierarchy is 
                    ℓ=0: $(sum( (all_splines₀[plot_idx].(xx_k) .- 𝒳₀_k[plot_idx].(xx_k)).^2 )/M),
                    ℓ=2: $(sum( (all_splines₂[plot_idx].(xx_k) .- 𝒳₂_k[plot_idx].(xx_k)).^2 )/M)\n")
            label="I = $(Nᵢ),"#M = $(M), lmax = $(ℓ_ν)" 
            plot!(p1,xx_k,abs.(𝒳₀_k[plot_idx].(xx_k)),label=label)
            plot!(p2,xx_k,abs.(𝒳₂_k[plot_idx].(xx_k)),label=false,c=p1.series_list[end][:linecolor])
            # plot!(p1,xx_k,𝒳₀_k[plot_idx].(xx_k),label=label)
            # plot!(p2,xx_k,𝒳₂_k[plot_idx].(xx_k),label=false,c=p1.series_list[end][:linecolor])
                    
        end
    end
    # ylims!(p1,-0.1,1.1)
    # ylims!(p2,-0.5,0.7)
    xlabel!(p2,L"x",xguidefontsize=18)
    ylabel!(p1,L"\mathcal{X}_{0}",xguidefontsize=18)
    ylabel!(p2,L"\mathcal{X}_{2}",xguidefontsize=18)
    l = @layout [a  ; b]
    title!(p1,"k = $(@sprintf("%.2f", ie_0_conf_late_c.ie.k/Mpcfac
    )), switch at $(@sprintf("%.1f", η_switch)) Mpc")
    p3 = plot(p1, p2, layout = l)
    savefig("../misc_plots/fft_debug/fft_experiments/mssv_k$(@sprintf("%.2f", ie_0_conf_late_c.ie.k/Mpcfac
            ))_q16log_planckmnu0p15_xiter_switch$(@sprintf("%.1f", η_switch)).pdf"
    )
# end
p3
η2x(η_switch/Mpcfac)