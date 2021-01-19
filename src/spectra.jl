
# # compute source functions and spectra
# struct SourceFunction{T, IT, AA<:AbstractArray}
#     S_grid::AA{T,2}
#     bessel::AA{IT,1}  # interpolator for bessel functions
# end


# OPTIMIZATION OPPORTUNITY
# should save u and du over the x_xgrid, it's an ODE option
# ℓᵧ is the Boltzmann hierarchy cutoff
function sourcefunction(par::AbstractCosmoParams{T}, bg, ih, k_grid,
        integrator::PerturbationIntegrator; ℓᵧ=8, reltol=1e-11) where T
    x_grid = bg.x_grid
    grid = zeros(T, length(x_grid), length(k_grid))
    @qthreads for (i_k, k) in enumerate(k_grid)
        hierarchy = Hierarchy(k, ℓᵧ, par, bg, ih)
        perturb = boltsolve(hierarchy, integrator; reltol=reltol)
        for (i_x, x) in enumerate(x_grid)
            u = perturb(x)
            du = similar(u)
            Bolt.basic_newtonian_hierarchy!(du, u, hierarchy, x)
            grid[i_x,i_k] = Bolt.source_function(du, u, hierarchy, x, hierarchy.par, integrator)
        end
    end
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
    bes = spline(bessel_xgrid, bessel_ygrid)
    return bes
end

function quadratic_k(kmin, kmax, nk)
    kmin, kmax, nk = assume_nondual(kmin), assume_nondual(kmax), assume_nondual(nk)
    return [kmin + (kmax - kmin) * (i/nk)^2 for i in 1:nk]
end

function Θl(k, s_itp, bes, par::AbstractCosmoParams{T}, bg) where {T}
    s = zero(T)
    xgrid = bg.x_grid
    for i in 1:length(xgrid)-1
        x = xgrid[i]
        sb = bes(k*(bg.η₀ - bg.η(x)))
        source = s_itp(x, k)
        s += sb * source * (xgrid[i+1] - xgrid[i])
    end
    return s
end

function cltt(ℓ, s_itp, kgrid, par::AbstractCosmoParams{T}, bg) where {T}
    bes = Bolt.bessel_interpolator(ℓ, kgrid[end] * bg.η₀)
    s = zero(T)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        th = Θl(k, s_itp, bes, par, bg)
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