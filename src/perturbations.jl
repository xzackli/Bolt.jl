
abstract type PerturbationIntegrator end
struct BasicNewtonian <: PerturbationIntegrator end

# purely for packing together everything needed to integrate a hierarchy at wavenumber k
struct Hierarchy{T<:Real, CP<:AbstractCosmoParams, BG<:AbstractBackground, IH<:AbstractIonizationHistory}
    k::T
    par::CP
    bg::BG
    ih::IH
end

# # convenience caller
# solve(k, par::AbstractCosmoParams, bg::AbstractBackground, ih::AbstractIonizationHistory,
#     integrator::PerturbationIntegrator) = solve(Hierarchy(k, par, bg, ih), integrator)


# basic newtonian gauge integrator
function solve(hierarchy::Hierarchy, integrator::BasicNewtonian, ode_alg; reltol=1e-10)
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = basic_newtonian_adiabatic_initial(hierarchy.par, hierarchy.bg, hierarchy.ih, xᵢ, hierarchy.k)
    prob = ODEProblem(basic_newtonian_hierarchy!, u₀, (xᵢ , 0.0), hierarchy)
    sol = solve(prob, ode_alg, reltol=reltol, dense=true)
    return sol
end

# specify a default integrator
solve(hierarchy::Hierarchy, integrator::BasicNewtonian; reltol=1e-10) = solve(
    hierarchy, integrator, Rodas5(); reltol=reltol)

function basic_newtonian_hierarchy!(du, u, hierarchy::Hierarchy, x)
    # compute cosmological quantities at time x, and do some unpacking
    k, par, bg, ih = hierarchy.k, hierarchy.par, hierarchy.bg, hierarchy.ih
    ℓᵧ, Ω_r, Ω_b, Ω_m, H₀² = par.ℓᵧ, par.Ω_r, par.Ω_b, par.Ω_m, bg.H₀^2
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.η(x), ih.τ′(x), ih.τ′′(x)

    # get array views of photon perturbations. this could be sugar'd up probably
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(view(du, 1:(ℓᵧ+1)), 0:ℓᵧ)
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(view(du, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]

    # preliminary variables
    a = x2a(x)
    R = 4Ω_r / (3Ω_b * a)

    # metric perturbations
    Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2])
    Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
        Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b + 4Ω_r * a^(-2) * Θ[0])

    # matter
    δ′ = k / ℋₓ * v - 3Φ′
    v′ = -v - k / ℋₓ * Ψ
    δ_b′ = k / ℋₓ * v_b - 3Φ′
        v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)

    # photons and polarized photons
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
    Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
    for ℓ in 2:(ℓᵧ-1)
        Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end

    Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
    for ℓ in 1:(ℓᵧ-1)
        Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
            (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
    end

    # photon boundary conditions: diffusion damping
    Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θ[ℓᵧ]
    Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * ηₓ) + τₓ′ * Θᵖ[ℓᵧ]

    du[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
    return nothing
end


function basic_newtonian_adiabatic_initial(par::AbstractCosmoParams{T,DT}, bg, ih, xᵢ, k) where {T,DT}
    ℓᵧ = par.ℓᵧ  # size of vector is determined by boltzmann hierarchy cutoff
    u = zeros(DT, 2ℓᵧ+7)
    ℋₓ, ℋₓ′, ηₓ, τₓ′, τₓ′′ = bg.ℋ(xᵢ), bg.ℋ′(xᵢ), bg.η(xᵢ), ih.τ′(xᵢ), ih.τ′′(xᵢ)
    Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ

    # metric and matter perturbations
    Φ = 1.0
    δ = 3Φ / 2
    δ_b = δ
    v = k / (2ℋₓ) * Φ
    v_b = v

    # photon hierarchy
    Θ[0] = Φ / 2
    Θ[1] = -k * Φ / (6ℋₓ)
    Θ[2] = -8k / (15ℋₓ * τₓ′) * Θ[1]
    Θᵖ[0] = (5/4) * Θ[2]
    Θᵖ[1] = -k / (4ℋₓ * τₓ′) * Θ[2]
    Θᵖ[2] = (1/4) * Θ[2]
    for ℓ in 3:ℓᵧ
        Θ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θ[ℓ-1]
        Θᵖ[ℓ] = -ℓ/(2ℓ+1) * k/(ℋₓ * τₓ′) * Θᵖ[ℓ-1]
    end

    u[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ, δ, v, δ_b, v_b  # pack back in
    return u
end


function basic_newtonian_source_function(du, u, hierarchy, x, par)
    # u = sol(x)
    # u′ = similar(u)
    # hierarchy!(u′, u, par, x)

    k, par, bg, ih = hierarchy.k, hierarchy.par, hierarchy.bg, hierarchy.ih
    ℓᵧ = par.ℓᵧ
    H₀² = bg.H₀^2
    ℋₓ, ℋₓ′, ℋₓ′′ = bg.ℋ(x), bg.ℋ′(x), bg.ℋ′′(x)
    τₓ, τₓ′, τₓ′′ = ih.τ(x), ih.τ′(x), ih.τ′′(x)
    g̃ₓ, g̃ₓ′, g̃ₓ′′ = ih.g̃(x), ih.g̃′(x), ih.g̃′′(x)

    # unpack variables from solution
    Θ = OffsetVector(u[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ = OffsetVector(u[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θ′ = OffsetVector(du[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Θᵖ′ = OffsetVector(du[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
    Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
    Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]
    Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]
    Φ′, δ′, v′, δ_b′, v_b′ = du[(2ℓᵧ+3):(2ℓᵧ+7)]

    a = x2a(x)
    Ψ = -Φ - 12H₀² / k^2 / a^2 * par.Ω_r * Θ[2]
    Ψ′ = -Φ′ - 12H₀² / k^2 / a^2 * par.Ω_r * (Θ′[2] - 2 * Θ[2])

    term1 =  g̃ₓ * (Θ[0] + Ψ + Π/4) + exp(-τₓ) * (Ψ′ - Φ′)
    term2 = (-1/k) * (ℋₓ′ * g̃ₓ * v_b + ℋₓ * g̃ₓ′ * v_b + ℋₓ * g̃ₓ * v_b′)
    Π′′ = 2k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * Θ[1] + Θ′[1]) + (3/10) * (τₓ′′ * Π + τₓ′ * Π′) -
        3k / (5ℋₓ) * (-ℋₓ′ / ℋₓ * (Θ[3] + Θᵖ[1] + Θᵖ[3]) + (Θ′[3] + Θᵖ′[1] + Θᵖ′[3]))
    term3 = (3/(4k^2)) * (
        (ℋₓ′^2 + ℋₓ * ℋₓ′′) * g̃ₓ * Π + 3 * ℋₓ * ℋₓ′ * (g̃ₓ′ * Π + g̃ₓ * Π′) +
        ℋₓ^2 * (g̃ₓ′′ * Π + 2g̃ₓ′ * Π′ + g̃ₓ * Π′′)
    )
    return term1 + term2 + term3
end



# ##

# par = CosmoParams()
# xᵢ = log(1e-8)
# xgridᵧ = collect(-20:0.005:0.0)
# H₀_ = H₀(par)
# zgrid = x2z.(xgridᵧ)
# η = Bolt.η_function(xgridᵧ, par)
# η₀ = η(0.0)
# kmin, kmax = 0.1H₀_, 1000H₀_
# nk = 100
# kgridᵧ = [kmin + (kmax - kmin) * (i/nk)^2 for i in 1:nk]

# ##
# using ThreadPools

# function generate_s_grid(par::AbstractCosmoParams{T,DT}, xgrid, kgrid) where {T,DT}
#     grid = zeros(DT, length(xgrid), length(kgrid))
#     @qthreads for k_i in eachindex(kgrid)
#         grid[:,k_i] .= source_x_grid(kgrid[k_i], xgrid, par)
#     end
#     return grid
# end

# @time s_kx_grid = generate_s_grid(par, xgridᵧ, kgridᵧ)

# ##
# using Interpolations
# # @time ss = source_x_grid(340H₀(par), xgrid, par)
# s_itp = LinearInterpolation((xgridᵧ, kgridᵧ), s_kx_grid, extrapolation_bc = Line())


# ℓ̂ = 100
# bessel_argmin = 0.0
# bessel_argmax = kmax * η₀
# Δg = bessel_argmax / 5000

# bessel_xgrid = bessel_argmin:Δg:bessel_argmax
# bessel_ygrid = [sphericalbesselj(ℓ̂, x) for x in bessel_xgrid]
# bes = LinearInterpolation((bessel_xgrid), bessel_ygrid, extrapolation_bc = Line())

# clf()
# k̂ = 340H₀(par)
# plot(xgridᵧ, [s_itp(x, k̂) * bes(k̂*(η₀ - η(x)))/1e-3 for x in xgridᵧ], "-", lw=0.5)

# ylabel(raw"Source function $\times$ bessel")
# xlabel(raw"$x$")
# xlim(-8, 0)
# ylim(-1, 3.5)
# # xlim(-8,-6)
# gcf()

# ##
# function Θl(k, s_itp, bes, xgrid, par::AbstractCosmoParams{T,DT}, η, η₀) where {T, DT}
#     s = zero(DT)
#     for i in 1:length(xgrid)-1
#         x = xgrid[i]
#         sb = bes(k*(η₀ - η(x)))::DT
#         source = s_itp(x, k)::DT
#         s += sb * source * (xgrid[i+1] - xgrid[i])
#     end
#     return s
# end

# ##
# @time Θl(340H₀_, s_itp, bes, xgridᵧ, par, η, η₀)
# ##

# clf()
# nk_dense = 5000
# dense_kgrid = [kmin + (kmax - kmin) * (i/nk_dense)^2 for i in 1:nk_dense]

# plot(dense_kgrid ./ H₀_,
#     [Θl(k, s_itp, bes, xgridᵧ, par, η, η₀)^2 / k / (1e-6 * H₀_^-1) for k in dense_kgrid],
#      "-")
# xlim(20,120)
# gcf()

# ##
# function Cl(ℓ, s_itp, xgrid, kgrid, par::AbstractCosmoParams{T,DT}, η, η₀) where {T,DT}
#     bessel_argmin = 0.0
#     bessel_argmax = kgrid[end] * η₀
#     Δg = bessel_argmax / 5000
#     bessel_xgrid = bessel_argmin:Δg:bessel_argmax
#     bessel_ygrid = [sphericalbesselj(ℓ, x) for x in bessel_xgrid]
#     bes = LinearInterpolation((bessel_xgrid), bessel_ygrid, extrapolation_bc = Line())

#     s = zero(DT)
#     for i in 1:length(kgrid)-1
#         k = kgrid[i]
#         dk = kgrid[i+1] - kgrid[i]
#         th = Θl(k, s_itp, bes, xgrid, par, η, η₀)::DT
#         s += th^2 * dk / k
#     end
#     return s
# end

# @time Cl(100, s_itp, xgridᵧ, dense_kgrid, par, η, η₀)


# ##
# function thCl(ells, s_itp, xgridᵧ, dense_kgrid, par, η, η₀)
#     cltt = zeros(length(ells))
#     @qthreads for (i,l) in enumerate(ells)
#         cltt[i] = Cl(l, s_itp, xgridᵧ, dense_kgrid, par, η, η₀)
#     end
#     return cltt
# end

# ells = 100:20:1200
# @time cltt = thCl(ells,  s_itp, xgridᵧ, dense_kgrid, par, η, η₀)

# ##
# clf()
# plt.plot(ells, cltt .* ells.^2)
# ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
# xlabel(raw"$\ell$")
# # yscale("log")
# gcf()

# ##
