#jms
# setup
using Revise
using Bolt
using Plots
using Interpolations
using DelimitedFiles
𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-25.0:0.1:0.0)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
x_grid = -25:0.1:0.0

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
# n_b(a, par) = par.Ω_b * ρ_crit(par) / (m_H * a^3)

# cross section
GF = 1.17e-5 * (1e-9)^2 #eV⁻¹, (1.17e-5GeV⁻¹)
σ_W(η) = GF^2 * T(η)^2

# decoupling time
Tdec = 0.8e6
adec = Tν0/Tdec
ηdec = a2η(adec)


# τ_ν′
τ_ν′(η) = n_e(η) * σ_W(η) * η2a(η) #/ bg.ℋ(log(η2a(η)))
ηmin = minimum(bg.η)
τ_ν′(ηmin),n_e(ηmin),η2a(ηmin),σ_W(ηmin)

# τ - do the integral using same setup as conformal time - checked this integrates to 1 as we ant
function τ_ν(η,η₀=bg.η(maximum(x_grid)))
  xq,wq = bg.quad_pts,bg.quad_wts
  logηmin,logηmax=log10(η),log10(η₀)
  #η = ∫ (aℋ)^-1 da/dx
  Iτ(y) = τ_ν′( xq2q(y,logηmin,logηmax) )  / dxdq(xq2q(y,logηmin,logηmax),logηmin,logηmax)
  return sum(Iτ.(xq).*wq)
end

# viz func
function g_ν(η,η₀ = bg.η(maximum(x_grid)))
  return τ_ν′(η) * exp(-τ_ν(η))
end

# check that viz func is a probability distribution in η
logηmint,logηmaxt=log10(minimum(bg.η)),log10(maximum(bg.η))
Ip(y) = g_ν( xq2q(y,logηmint,logηmaxt) )  / dxdq(xq2q(y,logηmint,logηmaxt),logηmint,logηmaxt)
sum(Ip.(bg.quad_pts).*bg.quad_wts)
#obviously this is not right but is "only" off by 30%

# viz func with factor of ℋ
plot(log10.(bg.η.(x_grid)),g_ν.(bg.η.(x_grid))./(bg.ℋ.(x_grid)),label=false)
vline!([log10(ηdec)],ls=:dot, label="naive decouple")
ylims!(0,1)
xlabel!(raw"$\log10 \eta$")
ylabel!(raw"$g_{\nu}(\eta) \mathcal{H}(\eta)^{-1}$")
ylims!(-0.001,1.1)
savefig("../nu_viz_func.png")

# add v_b term in 𝒩[1] proportional to τ_ν′

# source functions - do what I said and replace the dipole line with v_b
# make splines for τ,g -> add factor of ℋ to make tilde's to match current scfncs
spl_τ_ν = spline( [τ_ν(x2η(x))/  bg.ℋ(x) for x in x_grid], x_grid)
∂x_spl_τ_ν = spline_∂ₓ( spl_τ_ν, x_grid )
spl_g̃_ν = spline( [g_ν(x2η(x)) / bg.ℋ(x) for x in x_grid], x_grid)
∂x_spl_g̃_ν = spline_∂ₓ( spl_g̃_ν, x_grid )

function source_grid_ν(par::AbstractCosmoParams{T}, bg, ih, k_grid;
    ℓᵧ=10, reltol=1e-11,spl_τ_ν=spl_τ_ν,∂x_spl_τ_ν=∂x_spl_τ_ν,spl_g̃_ν=spl_g̃_ν,∂x_spl_g̃_ν=∂x_spl_g̃_ν) where T
    x_grid = bg.x_grid
    grid = zeros(T, length(x_grid), length(k_grid))
    for (i_k, k) in enumerate(k_grid)
        println("i_k = ", i_k)
        hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ)
        perturb = boltsolve(hierarchy; reltol=reltol)
        for (i_x, x) in enumerate(x_grid)
            u = perturb(x)  # this can be optimized away, save timesteps at the grid!
            du = similar(u)
            Bolt.hierarchy!(du, u, hierarchy, x;∂x_spl_τ_ν_x=∂x_spl_τ_ν(x)*bg.ℋ(x))
            grid[i_x,i_k] = source_function_ν(du, u, hierarchy, x,spl_τ_ν,spl_g̃_ν,∂x_spl_g̃_ν)
        end
    end
    # return grid
    itp = LinearInterpolation((x_grid, k_grid), grid, extrapolation_bc = Line())
    return itp
end

# CνB transfer
function 𝒩l(x_i, k, s_itp, bes, par::AbstractCosmoParams{T}, bg) where {T}
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

# c ells
function clnn(ℓ, s_itp, kgrid, par::AbstractCosmoParams{T}, bg) where {T}
    bes = Bolt.bessel_interpolator(ℓ, kgrid[end] * bg.η₀)
    x_i = findfirst(bg.x_grid .> -25)  #
    s = zero(T)
    for i in 1:length(kgrid)-1
        k = kgrid[i]
        dk = kgrid[i+1] - kgrid[i]
        N = 𝒩l(x_i, k, s_itp, bes, par, bg)
        s += N^2 * dk / k
    end
    println(ℓ)
    return s
end

function clnn(ℓ::Int, par::AbstractCosmoParams, bg, ih, sf)
    dense_kgrid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 5000)
    clnn(ℓ, sf, dense_kgrid, par, bg)
end

#---
#Actually run the Cl

𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-25.0:0.1:0.0)
# ih = IonizationHistory(Peebles(), 𝕡, bg)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
sf = source_grid_ν(𝕡, bg, ih, k_grid)
ℓs = 0:1:1200#500:520000
cells = [clnn(ℓ, 𝕡, bg, ih, sf) for ℓ in ℓs]

# save and plot
writedlm("../compare/fixRnu_coarse_vb2_bolt_cnb_cell.txt", (ℓs,cells))

plot(ℓs,ℓs.^2 .* cells,label=false)#label="vb-coarse",ls=:solid)
xlims!(0,3000)
xlabel!(raw"$\ell$")
ylabel!(raw"$\ell^{2} C_{\ell}$")
xgridtest=collect(-25.0:0.1:0.0)


#looking at sound speed and horizon

plot(xgridtest,log10.(ih.csb².(xgridtest).^(1/2)))
function rs_ν(η)
  xq,wq = bg.quad_pts,bg.quad_wts
  logηmin,logηmax=log10(x2η(minimum(bg.x_grid)))-50,log10(η)
  #η = ∫ (aℋ)^-1 da/dx
  Ics(y) = ih.csb²( log(η2a(xq2q(y,logηmin,logηmax))) )^(1/2)  / dxdq(xq2q(y,logηmin,logηmax),logηmin,logηmax)
  return sum(Ics.(xq).*wq)
end
function rs_rel(η)
  xq,wq = bg.quad_pts,bg.quad_wts
  logηmin,logηmax=log10(x2η(minimum(bg.x_grid)))-50,log10(η)
  #η = ∫ (aℋ)^-1 da/dx
  Ics(y) = ( 1/3 )^(1/2)  / dxdq(xq2q(y,logηmin,logηmax),logηmin,logηmax)
  return sum(Ics.(xq).*wq)
end


plot(xgridtest,log10.(rs_ν.(x2η.(xgridtest))))
plot!(xgridtest,log10.(x2η.(xgridtest)))
plot(xgridtest,log10(.75*π) .+ log10.(x2η.(xgridtest)) .- log10.(rs_ν.(x2η.(xgridtest))))
hline!([log10(200)])

log10(x2η(maximum(bg.x_grid))) - 24

τ_ν(x2η(maximum(bg.x_grid)))


rs_try= rs_ν(x2η(log(adec))) * kmsMpc_to_eV *3e5^2 /100  / 𝕡.h^2
rs_tryrec= rs_ν(x2η(log(1e-3))) * kmsMpc_to_eV *3e5^2 /100  / 𝕡.h^2
rs_tryrel= rs_rel(x2η(log(adec))) * kmsMpc_to_eV *3e5^2 /100  / 𝕡.h^2


ℓpapprox_peak_ν = 0.75*π*bg.η₀*bg.H₀*3e5/(100*rs_tryrel) #this is still too big...

plot(xgridtest,.5*log10.(ih.csb².(xgridtest)))
plot(xgridtest,.5*log10.(ones(length(xgridtest)).*1/3))
plot!(xgridtest,.5*log10.(ones(length(xgridtest)).*1 ./(3 .*(1 .+ exp.(xgridtest).* (3𝕡.Ω_b) ./ (4𝕡.Ω_r) ))))

plot!(xgridtest,.5log10.(3e5^2  *exp.(xgridtest) .* ih.csb².(xgridtest)))
