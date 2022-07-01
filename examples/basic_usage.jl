using Bolt, Plots

# Assign cosmological parameters
𝕡 = CosmoParams(Ω_c = 0.3) # set kwargs like so to change the default values
# Compute expansion history quantities
bg = Background(𝕡)
# Compute ionization history (via RECFAST)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)

# Matter power spectrum
kmin,kmax,nk = 10bg.H₀,5000bg.H₀,32
ks = log10_k(kmin,kmax,nk) # k grid
pL = [plin(k,𝕡,bg,ih) for k in ks]
p1 = plot(ks, vcat(pL...), xscale=:log10, yscale=:log10)


# CMB Cᵀᵀ(ℓ)
ℓmin,ℓmax,nℓ = 2,20,1200
ℓs = ℓmin:ℓmax:nℓ
kmin,kmax,nk = 0.1bg.H₀, 1000bg.H₀, 100
ks = quadratic_k(kmin,kmax,nk)
sf = source_grid(𝕡, bg, ih, ks, BasicNewtonian()) # set up LOS source function interpolator
Cᵀᵀ = cltt(ℓs, 𝕡, bg, ih, sf)
p2 = plot(ℓs, @.(ℓs^2*Cᵀᵀ))
