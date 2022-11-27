using Interpolations
using Bolt
using Plots

𝕡 = CosmoParams(Ω_c = 0.3) # set kwargs like so to change the default values
# Compute expansion history quantities
bg = Background(𝕡)
# Compute ionization history (via RECFAST)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)

# Matter power spectrum
kmin,kmax,nk = 10bg.H₀, 20000bg.H₀, 128
ks = log10_k(kmin,kmax,nk) # k grid
n_q,ℓᵧ,ℓ_ν,ℓ_mν,x=15,10,10,10,0
@time pL = [plin(k,𝕡,bg,ih,n_q,ℓᵧ,ℓ_ν,ℓ_mν,x) for k in ks]
p1 = plot(ks, vcat(pL...), xscale=:log10, yscale=:log10, label = "Linear")

InterpPmm = Interpolations.interpolate( log10.(vcat(pL...)), BSpline(Cubic(Line(OnGrid()))))

x = LinRange(log10(first(ks)), log10(last(ks)), length(ks))
InterpPmm = scale(InterpPmm, x)
InterpPmm = Interpolations.extrapolate(InterpPmm, Line())

params_halofit = Bolt.halofit_params(InterpPmm)

pNL_takahashi2012 =  [Bolt.takahashi2012(InterpPmm, params_halofit, 𝕡.Ω_b+𝕡.Ω_c, k) for k in ks]
pNL_takabird =  [Bolt.takabird(InterpPmm, params_halofit, 𝕡.Ω_b+𝕡.Ω_c, 𝕡.Ω_r, 𝕡.Σm_ν, 𝕡.h, k) for k in ks]
n_q,ℓᵧ,ℓ_ν,ℓ_mν,x=15,10,10,10,0
pL_check, pNL_check = Bolt.plin_and_nonlin(ks,𝕡,bg,ih,n_q,ℓᵧ,ℓ_ν,ℓ_mν,x)

plot!(p1, ks, vcat(pNL_takahashi2012...), label = "Takahashi2012")
plot!(p1, ks, vcat(pNL_takabird...), label = "Takabird")

println(pNL_check == vcat(pNL_takabird...))
