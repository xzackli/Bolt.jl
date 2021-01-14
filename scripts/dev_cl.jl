using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using OffsetArrays
using OrdinaryDiffEq
import Bolt: H₀
using ForwardDiff
using SpecialFunctions

# utility methods
δ_kron(i, j) = (i == j) ? 1 : 0
∂ₓ(f, x) = ForwardDiff.derivative(f, x)

function source_x_grid(k, xgrid, par)
    Xₑ = Bolt.saha_peebles_recombination(par)
    τ, τ′, τ′′ = Bolt.τ_functions(xgrid, Xₑ, par)
    g̃, g̃′ = Bolt.g̃_functions(τ, τ′, τ′′)
    g̃′′ = x -> ∂ₓ(g̃′, x)
    η = Bolt.η_function(xgrid, par)
    ℋ, ℋ′ = x -> Bolt.ℋ(x, par), x -> Bolt.ℋ′(x, par)


    TCA_condition(k, ℋₓ, τₓ′) = (abs(k / (ℋₓ * τₓ′)) < 0.1) & (abs(τₓ′) > 10.0)

    # NOTE: NO NEUTRINOS 𝒩
    function hierarchy!(du, u, p::AbstractCosmoParams, x)
        ℓᵧ, Ω_r, Ω_b, Ω_m = p.ℓᵧ, p.Ω_r, p.Ω_b, p.Ω_m
        H₀² = H₀(p)^2
        ℋₓ, ℋₓ′ = ℋ(x), ℋ′(x)
        τₓ′, τₓ′′ = τ′(x), τ′′(x)
        a = x2a(x)
        R = 4Ω_r / (3Ω_b * a)

        # get array views of photon perturbations
        Θ = OffsetVector(view(u, 1:(ℓᵧ+1)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
        Θ′ = OffsetVector(view(du, 1:(ℓᵧ+1)), 0:ℓᵧ)
        Θᵖ = OffsetVector(view(u, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)  # indicies 0 through ℓᵧ
        Θᵖ′ = OffsetVector(view(du, (ℓᵧ+2):(2ℓᵧ+2)), 0:ℓᵧ)
        Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]

        # metric perturbations
        Ψ = -Φ - 12H₀² / k^2 / a^2 * (Ω_r * Θ[2])
        Φ′ = Ψ - k^2 / (3ℋₓ^2) * Φ + H₀² / (2ℋₓ^2) * (
            Ω_m * a^(-1) * δ + Ω_b * a^(-1) * δ_b + 4Ω_r * a^(-2) * Θ[0])

        # matter
        δ′ = k / ℋₓ * v - 3Φ′
        v′ = -v - k / ℋₓ * Ψ
        δ_b′ = k / ℋₓ * v_b - 3Φ′

        Θ′[0] = -k / ℋₓ * Θ[1] - Φ′
        if TCA_condition(k, ℋₓ, τₓ′)
            Θ′[2] = 0.0  # could be solved for
            term_3Θ′_vb′ = (
                -((1-R)*τₓ′ + (1+R)*τₓ′′) * (3Θ[1] + v_b) - k * Ψ / ℋₓ
                + (1 - ℋₓ′ / ℋₓ) * (k / ℋₓ) * (-Θ[0] + 2Θ[2]) + k / ℋₓ * (-Θ′[0] + 2Θ′[2])
            ) / ((1+R)*τₓ′ + ℋₓ′ / ℋₓ - 1)
            v_b′ = (-v_b - k * Ψ / ℋₓ + R * (
                term_3Θ′_vb′ + k / ℋₓ * (-Θ[0] + 2Θ[2]) - k / ℋₓ * Ψ
            )) / (1+R)
            Θ′[1] = (term_3Θ′_vb′ - v_b′) / 3
        else
            v_b′ = -v_b - k / ℋₓ * Ψ + τₓ′ * R * (3Θ[1] + v_b)
            Θ′[1] = k / (3ℋₓ) * Θ[0] - 2k / (3ℋₓ) * Θ[2] + k / (3ℋₓ) * Ψ + τₓ′ * (Θ[1] + v_b/3)
        end

        # photons
        Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
        for ℓ in 2:(ℓᵧ-1)
            Θ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ-1] -
                (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θ[ℓ+1] + τₓ′ * (Θ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
        end
        # polarized photons
        Θᵖ′[0] = -k / ℋₓ * Θᵖ[1] + τₓ′ * (Θᵖ[0] - Π / 2)
        for ℓ in 1:(ℓᵧ-1)
            Θᵖ′[ℓ] = ℓ * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ-1] -
                (ℓ+1) * k / ((2ℓ+1) * ℋₓ) * Θᵖ[ℓ+1] + τₓ′ * (Θᵖ[ℓ] - Π * δ_kron(ℓ, 2) / 10)
        end
        # photon hierarchy boundary conditions
        Θ′[ℓᵧ] = k / ℋₓ * Θ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * η(x)) + τₓ′ * Θ[ℓᵧ]
        Θᵖ′[ℓᵧ] = k / ℋₓ * Θᵖ[ℓᵧ-1] - (ℓᵧ + 1) / (ℋₓ * η(x)) + τₓ′ * Θᵖ[ℓᵧ]

        du[(2ℓᵧ+3):(2ℓᵧ+7)] .= Φ′, δ′, v′, δ_b′, v_b′  # put non-photon perturbations back in
        return nothing
    end


    function adiabatic_initial_conditions(par::AbstractCosmoParams{T,DT}, xᵢ) where {T,DT}
        ℓᵧ = par.ℓᵧ
        u = zeros(DT, 2ℓᵧ+7)
        ℋₓ = Bolt.ℋ(xᵢ, par)
        τₓ′ = Bolt.τ′(xᵢ, Xₑ, par)
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

    xᵢ = first(xgrid)
    u₀ = adiabatic_initial_conditions(par, xᵢ)
    prob = ODEProblem(hierarchy!, u₀, (xᵢ , 0.0), par)
    sol = solve(prob, AutoTsit5(Rodas5()), reltol=1e-10)

    function source_function(sol, k, x, par)
        u = sol(x)
        u′ = similar(u)
        hierarchy!(u′, u, par, x)

        ℓᵧ = par.ℓᵧ
        H₀² = H₀(par)^2
        ℋₓ, ℋₓ′, ℋₓ′′ = ℋ(x), ℋ′(x), ∂ₓ(ℋ′, x)
        τₓ, τₓ′, τₓ′′ = τ(x), τ′(x), τ′′(x)
        g̃ₓ, g̃ₓ′, g̃ₓ′′ = g̃(x), g̃′(x), g̃′′(x)

        # unpack variables from solution
        Θ = OffsetVector(u[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
        Θᵖ = OffsetVector(u[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
        Θ′ = OffsetVector(u′[1:(ℓᵧ+1)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
        Θᵖ′ = OffsetVector(u′[(ℓᵧ+2):(2ℓᵧ+2)], 0:ℓᵧ)  # indicies 0 through ℓᵧ
        Π = Θ[2] + Θᵖ[2] + Θᵖ[0]
        Π′ = Θ′[2] + Θᵖ′[2] + Θᵖ′[0]
        Φ, δ, v, δ_b, v_b = u[(2ℓᵧ+3):(2ℓᵧ+7)]
        Φ′, δ′, v′, δ_b′, v_b′ = u′[(2ℓᵧ+3):(2ℓᵧ+7)]

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


    S̃(x) = source_function(sol, k, x, par)
    return S̃.(xgrid)
end

##

par = CosmoParams()
xᵢ = log(1e-8)
xgridᵧ = collect(-20:0.005:0.0)
H₀_ = H₀(par)
zgrid = x2z.(xgridᵧ)
η = Bolt.η_function(xgridᵧ, par)
η₀ = η(0.0)
kmin, kmax = 0.1H₀_, 1000H₀_
nk = 100
kgridᵧ = [kmin + (kmax - kmin) * (i/nk)^2 for i in 1:nk]

##
using ThreadPools

function generate_s_grid(par::AbstractCosmoParams{T,DT}, xgrid, kgrid) where {T,DT}
    grid = zeros(DT, length(xgrid), length(kgrid))
    @qthreads for k_i in eachindex(kgrid)
        grid[:,k_i] .= source_x_grid(kgrid[k_i], xgrid, par)
    end
    return grid
end

@time s_kx_grid = generate_s_grid(par, xgridᵧ, kgridᵧ)

##
using Interpolations
# @time ss = source_x_grid(340H₀(par), xgrid, par)
s_itp = LinearInterpolation((xgridᵧ, kgridᵧ), s_kx_grid, extrapolation_bc = Line())


ℓ̂ = 100
bessel_argmin = 0.0
bessel_argmax = kmax * η₀
Δg = bessel_argmax / 5000

bessel_xgrid = bessel_argmin:Δg:bessel_argmax
bessel_ygrid = [sphericalbesselj(ℓ̂, x) for x in bessel_xgrid]
bes = LinearInterpolation((bessel_xgrid), bessel_ygrid, extrapolation_bc = Line())

clf()
k̂ = 340H₀(par)
plot(xgridᵧ, [s_itp(x, k̂) * bes(k̂*(η₀ - η(x)))/1e-3 for x in xgridᵧ], "-", lw=0.5)

ylabel(raw"Source function $\times$ bessel")
xlabel(raw"$x$")
xlim(-8, 0)
ylim(-1, 3.5)
# xlim(-8,-6)
gcf()

##
function Θl(k, s_itp, bes, xgrid, par::AbstractCosmoParams{T,DT}, η, η₀) where {T, DT}
    s = zero(DT)
    for i in 1:length(xgrid)-1
        x = xgrid[i]
        sb = bes(k*(η₀ - η(x)))::DT
        source = s_itp(x, k)::DT
        s += sb * source * (xgrid[i+1] - xgrid[i])
    end
    return s
end

##
@time Θl(340H₀_, s_itp, bes, xgridᵧ, par, η, η₀)
##

clf()
nk_dense = 5000
dense_kgrid = [kmin + (kmax - kmin) * (i/nk_dense)^2 for i in 1:nk_dense]

plot(dense_kgrid ./ H₀_,
    [Θl(k, s_itp, bes, xgridᵧ, par, η, η₀)^2 / k / (1e-6 * H₀_^-1) for k in dense_kgrid],
     "-")
xlim(20,120)
gcf()

##
function Cl(ℓ, s_itp, xgrid, kgrid, par::AbstractCosmoParams{T,DT}, η, η₀) where {T,DT}
    bessel_argmin = 0.0
    bessel_argmax = kgrid[end] * η₀
    Δg = bessel_argmax / 5000
    bessel_xgrid = bessel_argmin:Δg:bessel_argmax
    bessel_ygrid = [sphericalbesselj(ℓ, x) for x in bessel_xgrid]
    bes = LinearInterpolation((bessel_xgrid), bessel_ygrid, extrapolation_bc = Line())

    s = zero(DT)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        th = Θl(k, s_itp, bes, xgrid, par, η, η₀)::DT
        s += th^2 * dk / k
    end
    return s
end

@time Cl(100, s_itp, xgridᵧ, dense_kgrid, par, η, η₀)


##
function thCl(ells, s_itp, xgridᵧ, dense_kgrid, par, η, η₀)
    cltt = zeros(length(ells))
    @qthreads for (i,l) in enumerate(ells)
        cltt[i] = Cl(l, s_itp, xgridᵧ, dense_kgrid, par, η, η₀)
    end
    return cltt
end

ells = 100:20:1200
@time cltt = thCl(ells,  s_itp, xgridᵧ, dense_kgrid, par, η, η₀)

##
clf()
plt.plot(ells, cltt .* ells.^2)
ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
# yscale("log")
gcf()

##
