using Bolt
using ForwardDiff
using PyPlot
using BenchmarkTools


# Cₗ function of baryon density
function clb(Ω_b::DT, ells) where DT
    par = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(par)
    ih = IonizationHistory(Peebles(), par, bg)
    k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 500)
    sf = sourcefunction(par, bg, ih, k_grid, BasicNewtonian())
    return cltt(ells, par, bg, ih, sf)
end

ells = 10:2:1200
f(Ω_b) = clb(Ω_b, ells)
@time cl = f(0.046)
@time ∂cl = ForwardDiff.derivative(f, 0.046)
##
Δ = 1e-3
@time finitediff_∂cl = (f(0.046 + Δ) .- f(0.046 - Δ)) ./ 2Δ
##

##
clf()
plt.plot(ells, cl .* ells.^2, "-", label=raw"$C_{\ell}$")
plt.plot(ells, ∂cl .* ells.^2 / 10, "-",
    label=raw"$\partial C_{\ell}/\partial\Omega_b / 10$")
# plt.plot(ells, finitediff_∂cl .* ells.^2 / 10, "--",
#     label=raw"$\Delta C_{\ell}/\Delta\Omega_b / 10$")

ylabel(raw"$\ell^2 C_{\ell}^{TT}$")
xlabel(raw"$\ell$")
legend()
ylim(-0.3, 0.5)
savefig("docs/assets/example_spectrum.png")
gcf()
