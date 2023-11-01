
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

function cltt(ℓ, s_itp, kgrid, par::AbstractCosmoParams{T}, bg) where {T}
    bes = Bolt.bessel_interpolator(ℓ, kgrid[end] * bg.η₀)
    x_i = findfirst(bg.x_grid .> -8)  # start integrating after recombination
    s = zero(T)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        th = Θl(x_i, k, s_itp, bes, par, bg)
        s += th^2 * dk / k
    end
    return s
end

function cltt(ℓ::Int, par::AbstractCosmoParams, bg, ih, sf)
    dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000)
    cltt(ℓ, sf, dense_kgrid, par, bg)
end

function cltt(ℓ⃗, par::AbstractCosmoParams, bg, ih, sf)
    dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000)
    return qmap(ℓ->cltt(ℓ, par, bg, ih, sf), ℓ⃗)
end


function plin(k, 𝕡::AbstractCosmoParams{T},bg,ih,ℓᵧ=15,x=0,reltol=1e-5) where T
    #copy code abvoe
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ) #shoddy quality test values
    perturb = boltsolve(hierarchy; reltol=reltol)
    results = perturb(x)
    #Also using the fact that a=1 at z=0
    δcN,δbN = results[2(ℓᵧ+1)+2,:],results[2(ℓᵧ+1)+4,:]
    vcN,vbN = results[2(ℓᵧ+1)+3,:],results[2(ℓᵧ+1)+5,:]
    #omegas to get weighted sum for total matter in background
    Ωm = 𝕡.Ω_c+𝕡.Ω_b
    #construct gauge-invariant versions of density perturbations
    δc = δcN - 3bg.ℋ(x)*vcN /k 
    δb = δbN - 3bg.ℋ(x)*vbN /k 
    δm = (𝕡.Ω_c*δc + 𝕡.Ω_b*δb) / Ωm
    As=𝕡.A
    Pprim = As*(k/0.05)^(𝕡.n-1) #pivot scale from Planck (in Mpc^-1)
    PL= (2π^2 / k^3)*δm^2 *Pprim
    return PL
end
