using Bolt, Plots, DelimitedFiles

#FIXME adjust the cl test accordingly

# comparison data
prefix="./test/data/camb_v_class_2022_default_camb_cmb_"
# reduce the number of ells since it is too many
ℓ_stride = 40
camb_Cᵀᵀ = readdlm(prefix*"unlens_cellTT.dat")[3:ℓ_stride:end,2] #start at 3 because that is first nonzero ℓ
camb_Cᵀᴱ = readdlm(prefix*"unlens_cellTE.dat")[3:ℓ_stride:end,2]
camb_Cᴱᴱ = readdlm(prefix*"unlens_cellEE.dat")[3:ℓ_stride:end,2]
camb_ℓs = readdlm(prefix*"unlens_cellTT.dat")[3:ℓ_stride:end,1]
ℓs = Int.(camb_ℓs); #FIXME silent error issues if Float64

# global ℓ range
# ℓmin,ℓmax,nℓ = 2,20,1200
# ℓs = ℓmin:ℓmax:nℓ

# subsample the ℓs because 2401 is too many
# (TODO check how camb gets so many ℓs??)

function clttteee(ℓs, 𝕡, bg, ih, sf,sf_P; Dℓ=true)
    Cᵀᵀ = cltt(ℓs, 𝕡, bg, ih, sf)
    Cᵀᴱ = clte(ℓs, 𝕡, bg, ih, sf,sf_P)
    Cᴱᴱ = clee(ℓs, 𝕡, bg, ih, sf_P)
    Tγ = (15/ π^2 *Bolt.ρ_crit(𝕡) * 𝕡.Ω_r)^(1/4)
    norm_fac = ( (Tγ * Bolt.Kelvin_natural_unit_conversion * 1e6)^2 / (2π) * 4^2)
    if Dℓ
        ℓsq = ℓs.^2 *norm_fac
        return ℓsq.*Cᵀᵀ, ℓsq.*Cᵀᴱ, ℓsq.*Cᴱᴱ
    else
        return Cᵀᵀ*norm_fac, Cᵀᴱ*norm_fac, Cᴱᴱ*norm_fac
    end
end


# test this function
function setup_sf(nx,nk)
    # setup
    𝕡 = CosmoParams()
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

function single_experiment(nx,nk)
    # setup
    𝕡, bg, ih, sf, sf_P = setup_sf(nx,nk)

    # compute Cℓs
    Cᵀᵀ, Cᵀᴱ, Cᴱᴱ = clttteee(ℓs, 𝕡, bg, ih, sf,sf_P)

    # compare to CAMB
    rᵀᵀ,rᵀᴱ,rᴱᴱ = Cᵀᵀ.-camb_Cᵀᵀ, Cᵀᴱ.-camb_Cᵀᴱ, Cᴱᴱ.-camb_Cᴱᴱ

    return rᵀᵀ,rᵀᴱ,rᴱᴱ
end

# CMB Cᵀᵀ(ℓ)
function plot_experiments(rrs,nxs,nks,percent=false,sv=["TT","TE","EE"])
    p = hline([0.0],color=:black,ls=:solid,label=false)
    for (i,rr) in enumerate(rrs)
        for s in sv
            if percent
                if s == "TT"
                    camb_C = camb_Cᵀᵀ
                elseif s == "TE"
                    camb_C = camb_Cᵀᴱ
                elseif s == "EE"
                    camb_C = camb_Cᴱᴱ
                else
                    error("Invalid spectrum")
                end
                label_use = s*" - nx = $(nxs[i]), nk = $(nks[i])" 
                plot!(ℓs, rr[1]./camb_C, label=label_use,ls=:solid,color=i)
            else
                label_use = s*" - nx = $(nxs[i]), nk = $(nks[i])" 
                plot!(ℓs, rr[1], label=label_use,ls=:solid,color=i)
            end
        end
    end
    xlabel!(p,raw"$\ell$")
    if percent
        ylabel!(p,raw"$\Delta C_{\ell}^{(\mathrm{CAMB})} / C_{\ell}^{(\mathrm{CAMB})}$")
        ylims!(-0.5,0.5)
    else
        ylabel!(p,raw"$\Delta C_{\ell}^{(\mathrm{CAMB})}$")
    end
    title!(p,"Scalar Cℓ convergence wrt. CAMB - percent = $percent")
    return p
end

nxx = [100,200,300,400]
nkk = [50,100,150,200]

# Structured runs
# the scheme is a list with [(nx_1,nk_1), (nx_1,nk_2),...,(nx_2,nk_1),...]
rrs = [single_experiment(nx,nk) 
        for nk in nkk 
         for nx in nxx]

nxx_2d = [nxx[i_k] for i_k in 1:length(nkk) for i_x in 1:length(nxx)]
nkk_2d = vcat([nkk for i_x in 1:length(nxx)]...)


# can take "TT","TE","EE" or any combination, but gets cluttered
p = plot_experiments(rrs,nxx_2d,nkk_2d,true,["TT"]) 
savefig(p,save_path*"rel_cltt_diffs.png")
p = plot_experiments(rrs,nxx_2d,nkk_2d,true,["TE"]) 
savefig(p,save_path*"rel_clte_diffs.png")
p = plot_experiments(rrs,nxx_2d,nkk_2d,true,["EE"]) 
savefig(p,save_path*"rel_clee_diffs.png")
# save the results
save_path = "../../plots/"
writedlm(save_path*"clttteee_diffs.dat",rrs)



