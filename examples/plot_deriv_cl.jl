using Bolt
using ForwardDiff
# using PyPlot
using Plots
using BenchmarkTools

# Cₗ as a function of baryon density
#Some noise at lowest few ells...9/17/21

function clb(Ω_b::DT, ells) where DT
    𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
    sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    return cltt(ells, 𝕡, bg, ih, sf)
end

ells = 10:20:1200
f(Ω_b) = clb(Ω_b, ells)
@time cl = f(0.046)
@time ∂cl = ForwardDiff.derivative(f, 0.046)  # you can just ForwardDiff the whole thing

##
# clf()
# plt.figure(figsize=(10,5))
plot(ells, cl .* ells.^2, label=raw"$C_{\ell}$",lw=0,marker=:circle)
plot!(ells, ∂cl .* ells.^2 / 10,
    label=raw"$\partial C_{\ell}/\partial\Omega_b / 10$",lw=0,marker=:circle)

# uncomment to see finite differences
# Δ = 1e-3
@time finitediff_∂cl = (f(0.046 + Δ) .- f(0.046 - Δ)) ./ 2Δ
plot!(ells, finitediff_∂cl .* ells.^2 / 10, ls=:dash,
    label=raw"$\Delta C_{\ell}/\Delta\Omega_b / 10$")

ylabel!(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel!(raw"$\ell$")
ylims!(-0.3, 0.5)
# savefig("docs/assets/example_spectrum.png")
