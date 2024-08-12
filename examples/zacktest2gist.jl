using Bolt, Plots, DelimitedFiles

#FIXME adjust the cl test accordingly

# comparison data
prefix="./test/data/camb_v_class_2022_default_camb_cmb_"
# reduce the number of ells since it is too many
ℓ_stride = 40
𝕡 = CosmoParams();
Tγ = (15/ π^2 *Bolt.ρ_crit(𝕡) * 𝕡.Ω_r)^(1/4);
norm_fac = (Tγ*Bolt.Kelvin_natural_unit_conversion * 1e6)^2;
norm_fac
# check this is the same cmb temp as camb generated these files...
camb_Cᵀᵀ = readdlm(prefix*"unlens_cellTT.dat")[3:ℓ_stride:end,2]/norm_fac; #start at 3 because that is first nonzero ℓ
camb_Cᵀᴱ = readdlm(prefix*"unlens_cellTE.dat")[3:ℓ_stride:end,2]/norm_fac;
camb_Cᴱᴱ = readdlm(prefix*"unlens_cellEE.dat")[3:ℓ_stride:end,2]/norm_fac;
camb_ℓs = readdlm(prefix*"unlens_cellTT.dat")[3:ℓ_stride:end,1];
ℓs = Int.(camb_ℓs); #FIXME silent error issues if Float64


function setup_sf(nx,nk)
    # setup
    # 𝕡 = CosmoParams()
    xmin,xmax = -20.0,0.0
    dx = (xmax-xmin)/nx
    bg = Background(𝕡; x_grid=xmin:dx:xmax)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    kmin,kmax = 0.1bg.H₀, 1000bg.H₀ #WARNING fd issues with these 𝕡-dep. values
    ks = quadratic_k(kmin,kmax,nk)
    sf = source_grid(𝕡, bg, ih, ks, BasicNewtonian())
    sf_P = source_grid_P(𝕡, bg, ih, ks, BasicNewtonian())
    return 𝕡, bg, ih, sf, sf_P
end

𝕡, bg, ih, sf, sf_P = setup_sf(400,200)

##

x_i_test = findfirst(bg.x_grid .> -8)  # start integrating after recombination
ks_test = quadratic_k(0.1bg.H₀,1000bg.H₀,200);

tΘl_arr = zeros(length(camb_ℓs),length(ks_test));
for (i,ℓ) in enumerate(camb_ℓs)
    bes = Bolt.bessel_interpolator(ℓ, maximum(ks_test) * bg.η₀)
    for (j,k) in enumerate(ks_test[1:end-1])
        kmid = (k + ks_test[j+1]) / 2
        tΘl_arr[i,j] = Bolt.Tl(x_i_test,kmid,sf,bes,bg)
    end
end

##
function t_cltt(ℓ_idx, kgrid) 
    s = zero(Float64)
    # bes = Bolt.bessel_interpolator(camb_ℓ[ℓ_idx], maximum(ks_test) * bg.η₀)
    for i in 1:length(kgrid)-1
        k = (kgrid[i] + kgrid[i+1])/2 #use midpoint
        dk = kgrid[i+1] - kgrid[i]
        # th = Bolt.Tl(x_i_test, k, sf, bes, bg) #
        th = tΘl_arr[ℓ_idx,i] 
        k_hMpc=k/(bg.H₀*2.99792e5/100) #This is messy...
        Pprim = 𝕡.A*(k_hMpc/0.05)^(𝕡.n-1)
        s += th^2 * Pprim * dk / k
    end
    return 4π*s
end


cltt_besp_arr = zeros(length(camb_ℓs))
cltt_arr = zeros(length(camb_ℓs))
for (i,ℓ) in enumerate(camb_ℓs)
    cltt_arr[i] = cltt(ℓ, sf, ks_test, 𝕡, bg)
    cltt_besp_arr[i] = t_cltt(i,ks_test)
end

plot(camb_ℓs,camb_Cᵀᵀ,label="CAMB")
scatter!(camb_ℓs,camb_ℓs.*(camb_ℓs.+1)/(2π) .* cltt_arr,label="Bolt")
scatter!(camb_ℓs,camb_ℓs.*(camb_ℓs.+1)/(2π) .* cltt_besp_arr,label="Bolt_besp", ms=2)
xlabel!(raw"$\ell$")
ylabel!(raw"$\ell(\ell+1)C_{\ell}^{TT}/2\pi$")