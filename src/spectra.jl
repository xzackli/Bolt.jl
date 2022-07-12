
# OPTIMIZATION OPPORTUNITY
# should save u and du over the x_xgrid, it's an ODE option
# ℓᵧ is the Boltzmann hierarchy cutoff
function source_grid(par::AbstractCosmoParams{T}, bg, ih, k_grid,
        integrator::PerturbationIntegrator; ℓᵧ=8, reltol=1e-11) where T
    x_grid = bg.x_grid
    grid = zeros(T, length(x_grid), length(k_grid))
    @qthreads for (i_k, k) in enumerate(k_grid)
        hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ)
        perturb = boltsolve(hierarchy; reltol=reltol)
        for (i_x, x) in enumerate(x_grid)
            u = perturb(x)  # this can be optimized away, save timesteps at the grid!
            du = similar(u)
            Bolt.hierarchy!(du, u, hierarchy, x)
            grid[i_x,i_k] = Bolt.source_function(du, u, hierarchy, x)
        end
    end
    # return grid
    itp = LinearInterpolation((x_grid, k_grid), grid, extrapolation_bc = Line())
    return itp
end

function source_grid_P(par::AbstractCosmoParams{T}, bg, ih, k_grid,
    integrator::PerturbationIntegrator; ℓᵧ=8, reltol=1e-11) where T
    x_grid = bg.x_grid
    grid = zeros(T, length(x_grid), length(k_grid))
    @qthreads for (i_k, k) in enumerate(k_grid)
        hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ)
        perturb = boltsolve(hierarchy; reltol=reltol)
        for (i_x, x) in enumerate(x_grid)
            u = perturb(x)  # this can be optimized away, save timesteps at the grid!
            du = similar(u)
            Bolt.hierarchy!(du, u, hierarchy, x)
            grid[i_x,i_k] = Bolt.source_function_P(du, u, hierarchy, x)
        end
    end
    # return grid
    itp = LinearInterpolation((x_grid, k_grid), grid, extrapolation_bc = Line())
    return itp
end

# we make the assumption that shifting the coordinates upon which we integrate
# does not affect our result. that is, we choose coordinates where the integral converges
assume_nondual(x::ForwardDiff.Dual) = ForwardDiff.value(x)
assume_nondual(x::Real) = x

function bessel_interpolator(ℓ, kmax_η₀)
    bessel_argmin = 0.0

    bessel_argmax = assume_nondual(kmax_η₀)
    Δg = bessel_argmax / 5000
    bessel_xgrid = bessel_argmin:Δg:bessel_argmax
    bessel_ygrid = [sphericalbesselj(ℓ, x) for x in bessel_xgrid]
    bes = spline(bessel_ygrid, bessel_xgrid)
    return bes
end

function quadratic_k(kmin::T, kmax::T, nk) where T
    kmin, kmax, nk = assume_nondual(kmin), assume_nondual(kmax), assume_nondual(nk)
    return T[kmin + (kmax - kmin) * (i/nk)^2 for i in 1:nk]
end

function log10_k(kmin::T, kmax::T, nk) where T
    kmin, kmax, nk = assume_nondual(kmin), assume_nondual(kmax), assume_nondual(nk)
    return T[10 ^(log10(kmin) + (log10(kmax/kmin))*(i-1)/(nk-1))  for i in 1:nk]
end

function Θl(x_i, k, s_itp, bes, par::AbstractCosmoParams{T}, bg) where {T}
    s = zero(T)
    xgrid = bg.x_grid
    for i in x_i:length(xgrid)-1
        x = xgrid[i]
        sb = bes(k*(bg.η₀ - bg.η(x)))
        source = s_itp(x, k)
        s += sb * source * (xgrid[i+1] - xgrid[i])
    end
    return s
end

# For polzn we have moved the factor of (bg.η₀ - bg.η(x))^-2 from S_p to here
function Pl(x_i, k, s_itp, bes, par::AbstractCosmoParams{T}, bg) where {T}
    s = zero(T)
    xgrid = bg.x_grid
    for i in x_i:length(xgrid)-1
        x = xgrid[i]
        sb = bes(k*(bg.η₀ - bg.η(x)))
        source = s_itp(x, k)
        s += sb * source * (xgrid[i+1] - xgrid[i])
    end
    return s
end

function cltt(ℓ, s_itp, kgrid, par::AbstractCosmoParams{T}, bg) where {T}
    bes = Bolt.bessel_interpolator(ℓ, kgrid[end] * bg.η₀)
    x_i = findfirst(bg.x_grid .> -8)  # start integrating after recombination
    s = zero(T)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        th = Θl(x_i, k, s_itp, bes, par, bg)
        k_hMpc=k/(bg.H₀*2.99792e5/100) #This is messy...
        Pprim = par.A*(k_hMpc/0.05)^(par.n-1)
        s += th^2 * Pprim * dk / k
    end
    return s
end

#jms 6/7/22 UNTESTED
function clte(ℓ, s_itp_t, s_itp_e, kgrid, par::AbstractCosmoParams{T}, bg) where {T}
    bes = Bolt.bessel_interpolator(ℓ, kgrid[end] * bg.η₀)
    x_i = findfirst(bg.x_grid .> -8)  # start integrating after recombination
    s = zero(T)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        th = Θl(x_i, k, s_itp_t, bes, par, bg)
        ep = Pl(x_i, k, s_itp_e, bes, par, bg) * ℓð
        k_hMpc=k/(bg.H₀*2.99792e5/100) #This is messy...
        Pprim = par.A*(k_hMpc/0.05)^(par.n-1)
        s += th * ep * Pprim * dk / k
    end
    return s
end

function clee(ℓ, s_itp_p, kgrid, par::AbstractCosmoParams{T}, bg) where {T}
    bes = Bolt.bessel_interpolator(ℓ, kgrid[end] * bg.η₀)
    ℓð = sqrt((ℓ+2)*(ℓ+1)*ℓ*(ℓ-1)) #spin factor
    x_i = findfirst(bg.x_grid .> -8)  # start integrating after recombination
    s = zero(T)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        ep = Pl(x_i, k, s_itp_p, bes, par, bg) * ℓð
        k_hMpc=k/(bg.H₀*2.99792e5/100) #This is messy...
        Pprim = par.A*(k_hMpc/0.05)^(par.n-1)
        s += ep^2 * Pprim *  dk / k
    end
    return s
end

function cltt(ℓ::Int, par::AbstractCosmoParams, bg, ih, sf)
    dense_kgrid = quadratic_k(0.01bg.H₀, 1000bg.H₀, 5000)
    cltt(ℓ, sf, dense_kgrid, par, bg)
end

function clte(ℓ::Int, par::AbstractCosmoParams, bg, ih, sf,sf_P)
    dense_kgrid = quadratic_k(0.01bg.H₀, 1000bg.H₀, 5000)
    clte(ℓ, sf, sf_P, dense_kgrid, par, bg)
end

function clee(ℓ::Int, par::AbstractCosmoParams, bg, ih,  sf_P)
    dense_kgrid = quadratic_k(0.01bg.H₀, 1000bg.H₀, 5000)
    clee(ℓ, sf_P, dense_kgrid, par, bg)
end

function cltt(ℓ⃗, par::AbstractCosmoParams, bg, ih, sf)
    # dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000) #what is the purpose of this?
    return qmap(ℓ->cltt(ℓ, par, bg, ih, sf), ℓ⃗)
end

function clte(ℓ⃗, par::AbstractCosmoParams, bg, ih, sf,sf_P)
    # dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000)
    return qmap(ℓ->clte(ℓ, par, bg, ih, sf,sf_P), ℓ⃗)
end

function clee(ℓ⃗, par::AbstractCosmoParams, bg, ih, sf_P)
    # dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000)
    return qmap(ℓ->clee(ℓ, par, bg, ih, sf_P), ℓ⃗)
end


function plin(k, 𝕡::AbstractCosmoParams{T},bg,ih,
              n_q=15,ℓᵧ=500,ℓ_ν=500,ℓ_mν=20,x=0) where T
    #copy code abvoe
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ,ℓ_ν,ℓ_mν,n_q) #shoddy quality test values
    perturb = boltsolve(hierarchy; reltol=1e-5)
    results = perturb(x)
    ℳρ,_ = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q],
                            bg,exp(x),𝕡)./ bg.ρ₀ℳ(x)
    #Below assumes negligible neutrino pressure for the normalization (fine at z=0)
    ℳθ = k*θ(results[2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+2n_q],
                     bg,exp(x),𝕡)./ bg.ρ₀ℳ(x)
    #Also using the fact that a=1 at z=0
    δcN,δbN = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h
    vcN,vbN = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+3,:],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5,:]* 𝕡.h
    ℳρN,ℳθN = ℳρ,ℳθ
    vmνN = -ℳθN ./ k
    #omegas to get weighted sum for total matter in background
    Tγ = (15/ π^2 *bg.ρ_crit *𝕡.Ω_r)^(1/4)
    ζ = 1.2020569
    νfac = (90 * ζ /(11 * π^4)) * (𝕡.Ω_r * 𝕡.h^2 / Tγ) *((𝕡.N_ν/3)^(3/4))
    #^the factor that goes into nr approx to neutrino energy density, plus equal sharing ΔN_eff factor for single massive neutrino
    Ω_ν = 𝕡.Σm_ν*νfac/𝕡.h^2
    Ωm = 𝕡.Ω_m+𝕡.Ω_b+Ω_ν

    #construct gauge-invariant versions of density perturbations
    δc = δcN - 3bg.ℋ(x)*vcN ./k
    δb = δbN - 3bg.ℋ(x)*vbN ./k
    #assume neutrinos fully non-relativistic and can be described by fluid (ok at z=0)
    δmν = ℳρN - 3bg.ℋ(x)*vmνN ./k
    δm = (𝕡.Ω_m*δc .+ 𝕡.Ω_b*δb .+ Ω_ν*δmν) ./ Ωm
    As=𝕡.A#1e-10*exp(3.043)
    k_hMpc=k/(bg.H₀*3e5/100)
    Pprim = As*(k_hMpc./0.05).^(𝕡.n-1)
    PL= (2π^2 ./ k_hMpc.^3).*(δm*𝕡.h).^2 .*Pprim
    return PL
end
