#jms
# setup
using Revise
using Bolt
using Plots
using Interpolations
𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-25.0:0.01:0.0)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
x_grid = -25:0.01:0.0

# pieces of g_ν

# conf time → a - probably not the most accurate...
x2η = spline( [bg.η(x) for x in x_grid], x_grid)
function η2a(η,dx=0.01)
  x_fine = minimum(x_grid):dx:maximum(x_grid)
  return exp(x_fine[argmin(abs.(x2η.(x_fine).-η))])
end
a2η(a) = bg.η(log(a))

# electron number density
const ζ = 1.2020569 #Riemann ζ(3) for phase space integrals
Tγ0 = (15/ π^2 *bg.ρ_crit *𝕡.Ω_r)^(1/4) #CMB temp today in eV
Tν0 = Tγ0 * (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3)
T(η) = Tν0 * η2a(η)^(-1) #ν temp going backwards from today's neutrino temp
n_e_rel(η) = ζ / (π^2) * 2 * 3/4 * T(η)^3
mp=938e6 #proton mass eV/c^2 (938 Mev)
n_b(η) = bg.ρ_crit*𝕡.Ω_b*η2a(η)^(-3) / mp
η_ani = a2η(Tν0/511e3) #around electron mass of 511keV
n_e(η) = (η<η_ani) ? n_e_rel(η) : n_b(η)
n_b(a, par) = par.Ω_b * ρ_crit(par) / (m_H * a^3)

# cross section
GF = 1.17e-5 * (1e-9)^2 #eV⁻¹, (1.17e-5GeV⁻¹)
σ_W(η) = GF^2 * T(η)^2

# decoupling time
Tdec = 0.8e6
adec = Tν0/Tdec
ηdec = a2η(adec)

# diagnostic plots
# plot(log10.(exp.(x_grid)),log10.(T.(bg.η(x_grid))))
# hline!([6])
# vline!([log10(2e-10)])

# τ_ν′
τ_ν′(η) = n_e(η) * σ_W(η) * η2a(η)
ηmin = minimum(bg.η)
τ_ν′(ηmin),n_e(ηmin),η2a(ηmin),σ_W(ηmin)

# diagnostic plots
# plot(log10.(bg.η.(x_grid)),log10.(τ_ν′.(bg.η.(x_grid))),ls=:dot)
# plot(log10.(bg.η.(x_grid)),log10.(n_e.(bg.η.(x_grid))))
# plot!(log10.(bg.η.(x_grid)),log10.(σ_W.(bg.η.(x_grid))))
# plot!(log10.(bg.η.(x_grid)),log10.(η2a.(bg.η.(x_grid))))
# vline!([log10(ηdec)],color=:black,ls=:dash)
# vline!([log10(η_ani)],color=:gray,ls=:dot)

# τ - do the integral using same setup as conformal time

function τ_ν(η,η₀=bg.η(maximum(x_grid)))
  xq,wq = bg.quad_pts,bg.quad_wts
  logηmin,logηmax=log10(η),log10(η₀)
  #η = ∫ (aℋ)^-1 da/dx
  Iτ(y) = τ_ν′( xq2q(y,logηmin,logηmax) )  / dxdq(xq2q(y,logηmin,logηmax),logηmin,logηmax)
  return sum(Iτ.(xq).*wq)
end

# viz func
function g_ν(η,η₀ = bg.η(maximum(x_grid)))
  return τ_ν′(η) / (bg.ℋ(log(η2a(η)))) * exp(-τ_ν(η))
end


# diagnostic plots

# plot(log10.(bg.η.(x_grid)),log10.(τ_ν.(bg.η.(x_grid))))
# plot(log10.(bg.η.(x_grid)),log10.(exp.(-τ_ν.(bg.η.(x_grid)))))
# plot!(log10.(bg.η.(x_grid)),log10.(τ_ν′.(bg.η.(x_grid))))
# plot!(log10.(bg.η.(x_grid)),log10.(τ_ν′.(bg.η.(x_grid))./(bg.ℋ.(x_grid))))
# plot!(log10.(bg.η.(x_grid)),log10.(g_ν.(bg.η.(x_grid))))

# viz func with factor of ℋ
plot(log10.(bg.η.(x_grid)),g_ν.(bg.η.(x_grid)),label=false)
vline!([log10(ηdec)],ls=:dot, label="naive decouple")
ylims!(0,1)
xlabel!(raw"$\log10 \eta$")
ylabel!(raw"$g_{\nu}(\eta) \mathcal{H}(\eta)^{-1}$")
ylims!(-0.001,1.1)
savefig("../nu_viz_func.png")

# v_b term

# source functions
