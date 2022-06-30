using Bolt
using ForwardDiff
using Plots
using BenchmarkTools

##
# P(k) as a function of Ω_c

function pkc(Ω_c::DT,k_grid) where DT
    𝕡 = CosmoParams{DT}(Ω_c=Ω_c)
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=1)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    # sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    nk=length(k_grid)
    a=zeros(DT,nk)
    for i = 1:nk
            #   pl = plin(k_grid[i],𝕡,bg,ih)
            #   println(i)
              a[i] = plin(k_grid[i],𝕡,bg,ih)[1] #not sure why this is returning a vector and not a T...
          end
   return a
end

#clean this up by specifying k in h/Mpc
𝕡 = CosmoParams(Ω_c=0.224)
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=15)
kmin,kmax= 0.1bg.H₀*100,5000bg.H₀
k_grid = log10_k(kmin,kmax,33)
k_grid_hMpc = k_grid/(bg.H₀*3e5/100)
fc(Ω_c) = pkc(Ω_c,k_grid)
@btime pk = fc(0.224)
@btime ∂pk = ForwardDiff.derivative(fc, 0.224)  # you can just ForwardDiff the whole thing
pk
##
∂pk.>0
pos_∂ = ∂pk[∂pk.>0]
neg_∂ = ∂pk[∂pk.<0]
pos_k = k_grid_hMpc[∂pk.>0]
neg_k = k_grid_hMpc[∂pk.<0]

plot(log10.(k_grid_hMpc), log10.(pk), label=raw"$P_{L}(k)$",lw=0,marker=:circle)
plot!(log10.(pos_k), log10.(pos_∂),
    label=raw"$\log ~ \partial P_{L}(k)/\partial \Omega_{c}$",
    lw=0,marker=:circle)
plot!(log10.(neg_k), log10.(-neg_∂),
    label=raw"$\log ~ -\partial P_{L}(k)/\partial \Omega_{c}$",
    lw=0,marker=:diamond,color=palette(:default)[2])
# uncomment to see finite differences
Δ = 1e-3
@time finitediff_∂pk = (fc(0.224 + Δ) .- fc(0.224 - Δ)) ./ 2Δ
plot!(log10.(k_grid_hMpc), -finitediff_∂pk ,  ls=:dash,
    label=raw"$ -\log ~ \Delta P_{L}(k)/ \Delta \Omega_{c}$")
ylabel!(raw"$\log ~ P_{L}(k) \ [\mathrm{Mpc}/h]^{3}$")
xlabel!(raw"$\log ~  k \ [h/\mathrm{Mpc}]$")
savefig("docs/assets/example_linear_power_c.png")

##
# P(k) as a function of neutrino mass

function pkm(Σmν::DT,k_grid) where DT
    𝕡 = CosmoParams{DT}(Σm_ν=Σmν)
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    # sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    nk=length(k_grid)
    a=zeros(DT,nk)
    for i = 1:nk
              pl = plin(k_grid[i],𝕡,bg,ih)
              println(i)
              a[i] = plin(k_grid[i],𝕡,bg,ih)[1] #not sure why this is returning a vector and not a T...
          end
   return a
end

#clean this up by specifying k in h/Mpc
𝕡 = CosmoParams(Σm_ν=0.1)
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
kmin,kmax= 0.1bg.H₀*100,5000bg.H₀
k_grid = log10_k(kmin,kmax,33)
k_grid_hMpc = k_grid/(bg.H₀*3e5/100)
f(Σmν) = pkm(Σmν,k_grid)
@time pk = f(0.1)
@time ∂pk = ForwardDiff.derivative(f, 0.1)  # you can just ForwardDiff the whole thing


##
# clf()
# plt.figure(figsize=(10,5))
plot(log10.(k_grid_hMpc), log10.(pk), label=raw"$P_{L}(k)$",lw=0,marker=:circle)
plot!(log10.(k_grid_hMpc), log10.(-∂pk),
    label=raw"$ -\log ~ \partial P_{L}(k)/\partial \Sigma m_{\nu}$",
    lw=0,marker=:circle)

# uncomment to see finite differences
# Δ = 1e-3
# @time finitediff_∂pk = (f(0.1 + Δ) .- f(0.1 - Δ)) ./ 2Δ
# plot!(log10.(k_grid_hMpc), log10.(-finitediff_∂pk) ,  ls=:dash,
#     label=raw"$ -\log ~ \Delta P_{L}(k)/ \Delta \Sigma m_{\nu}$")

ylabel!(raw"$\log ~ P_{L}(k) \ [\mathrm{Mpc}/h]^{3}$")
xlabel!(raw"$\log ~  k \ [h/\mathrm{Mpc}]$")
# savefig("docs/assets/example_linear_power_n.png")
