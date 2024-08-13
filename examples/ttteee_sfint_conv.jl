using Bolt, Plots, DelimitedFiles, Revise

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
camb_Cᵀᵀ = readdlm(prefix*"unlens_cellTT.dat")[3:ℓ_stride:end,2]/norm_fac #start at 3 because that is first nonzero ℓ
camb_Cᵀᴱ = readdlm(prefix*"unlens_cellTE.dat")[3:ℓ_stride:end,2]/norm_fac
camb_Cᴱᴱ = readdlm(prefix*"unlens_cellEE.dat")[3:ℓ_stride:end,2]/norm_fac
camb_ℓs = readdlm(prefix*"unlens_cellTT.dat")[3:ℓ_stride:end,1]
ℓs = Int.(camb_ℓs); #FIXME silent error issues if Float64
# subsample the ℓs because 2401 is too many
# (TODO check how camb gets so many ℓs??)


function clttteee(ℓs, 𝕡, bg, ih, sf,sf_P; Dℓ=true)
    Cᵀᵀ = cltt(ℓs, 𝕡, bg, ih, sf)
    Cᵀᴱ = clte(ℓs, 𝕡, bg, ih, sf,sf_P)
    Cᴱᴱ = clee(ℓs, 𝕡, bg, ih, sf_P)
    if Dℓ
        ℓsq = ℓs.*(ℓs.+1) ./ 2π # Dℓ = ℓ(ℓ+1)/(2π) C_ℓ
        return  ℓsq.*Cᵀᵀ, ℓsq.*Cᵀᴱ, ℓsq.*Cᴱᴱ
    else
        return Cᵀᵀ*norm_fac, Cᵀᴱ*norm_fac, Cᴱᴱ*norm_fac
    end
end


# test this function
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

function single_experiment(nx,nk)
    # setup
    𝕡, bg, ih, sf, sf_P = setup_sf(nx,nk)

    # compute Cℓs
    Cᵀᵀ, Cᵀᴱ, Cᴱᴱ = clttteee(ℓs, 𝕡, bg, ih, sf,sf_P)

    # compare to CAMB
    rᵀᵀ,rᵀᴱ,rᴱᴱ = Cᵀᵀ.-camb_Cᵀᵀ, Cᵀᴱ.-camb_Cᵀᴱ, Cᴱᴱ.-camb_Cᴱᴱ

    return rᵀᵀ,rᵀᴱ,rᴱᴱ,Cᵀᵀ,Cᵀᴱ,Cᴱᴱ,sf,sf_P
end

# CMB Cᵀᵀ(ℓ)
function plot_experiments(rrs,nxs,nks,nxs_1d,percent=false,sv=["TT","TE","EE"])
    p = hline([0.0],color=:black,ls=:solid,label=false,legend_columns=2,legend=percent ? :bottomright : :topright)
    for (i,rr) in enumerate(rrs)
        for s in sv
            if s == "TT"
                camb_C = camb_Cᵀᵀ
                r_idx = 1
            elseif s == "TE"
                camb_C = camb_Cᵀᴱ
                r_idx = 2
            elseif s == "EE"
                camb_C = camb_Cᴱᴱ
                r_idx = 3
            else
                error("Invalid spectrum")
            end
            if percent
                label_use = s*" - nx = $(nxs[i]), nk = $(nks[i])" 
                plot!(ℓs, rr[r_idx]./camb_C, label=label_use,ls=:solid,color=(i-1)%length(nxs_1d)+1,lw=(1+i÷length(nxs_1d)))
                println(i%length(nxs_1d))
            else
                label_use = s*" - nx = $(nxs[i]), nk = $(nks[i])" 
                plot!(ℓs, rr[r_idx], label=label_use,ls=:solid,color=(i-1)%length(nxs_1d)+1,lw=(1+i÷length(nxs_1d)))
                # plot!(ℓs,camb_C,ls=:solid,color=:black,label= i==1 ? "CAMB" : false)
            end
        end
    end
    xlabel!(p,raw"$\ell$")
    if percent
        ylabel!(p,raw"$\Delta C_{\ell} / C_{\ell}^{(\mathrm{CAMB})}$")
        # ylims!(-0.5,0.5)
    else
        ylabel!(p,raw"$\Delta C_{\ell}^{(\mathrm{CAMB})}$")
    end
    title!(p,"Scalar Cℓ convergence wrt. CAMB - percent = $percent")
    return p
end

nxx = [200,300,400]
nkk = [100,150,200]


all_test = single_experiment(400,200) 
bg_test = Background(𝕡; x_grid=-20.0:(0.0+20.0)/400:0.0);
ks_test = quadratic_k(0.1bg_test.H₀,1000bg_test.H₀,200);

# transfer fucntions
camb_T_ℓ_k = readdlm("./test/data/camb_v_class_2022_default_camb_cmb_transfers_T_ell_k.dat")
camb_E_ℓ_k = readdlm("./test/data/camb_v_class_2022_default_camb_cmb_transfers_E_ell_k.dat")
camb_ℓ = readdlm("./test/data/camb_v_class_2022_default_camb_cmb_transfers_ell.dat")
camb_k = readdlm("./test/data/camb_v_class_2022_default_camb_cmb_transfers_k.dat")

test_xi = findfirst(bg_test.x_grid .> -8)
function tΘl(x_i, k, s_itp, bes, bg)
    s = zero(Float64)
    xgrid = bg.x_grid
    for i in x_i:length(xgrid)-1
        x = xgrid[i]
        sb = bes(k*(bg.η₀ - bg.η(x)))
        source = s_itp(x, k)
        s += sb * source * (xgrid[i+1] - xgrid[i])
    end
    return s
end

tΘl_arr = zeros(length(camb_ℓ),length(ks_test))
for (i,ℓ) in enumerate(camb_ℓ)
    println(i)
    for (j,k) in enumerate(ks_test)
        bes = Bolt.bessel_interpolator(ℓ, k * bg_test.η₀)
        tΘl_arr[i,j] = tΘl(test_xi,k,all_test[7],bes,bg_test)
    end
end

#ℓ = 2
plot(camb_k, camb_T_ℓ_k[1,:],xscale=:log10,label="camb")
scatter!(ks_test, tΘl_arr[1,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$\Theta_{\ell=2}(k)$")
savefig("../../misc_plots/Θl_ℓ=2.png")

#ℓ = 40
plot(camb_k, camb_T_ℓ_k[26,:],xscale=:log10,label="camb")
scatter!(ks_test, tΘl_arr[26,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$\Theta_{\ell=40}(k)$")
savefig("../../misc_plots/Θl_ℓ=40.png")

#ℓ = 200
plot(camb_k, camb_T_ℓ_k[41,:],xscale=:log10,label="camb")
scatter!(ks_test, tΘl_arr[41,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$\Theta_{\ell=200}(k)$")
savefig("../../misc_plots/Θl_ℓ=200.png")

#ℓ = 2000
plot(camb_k, camb_T_ℓ_k[end-10,:],xscale=:log10,label="camb")
scatter!(ks_test, tΘl_arr[end-10,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$\Theta_{\ell=2000}(k)$")
savefig("../../misc_plots/Θl_ℓ=2000.png")

# polzn
function tEl(x_i, k, s_itp, bes, bg)
    # function Pl(x_i, k, s_itp, bes, bg) where {T}
        s = zero(Float64)
        xgrid = bg.x_grid
        for i in x_i:length(xgrid)-1
            x = xgrid[i]
            sb = bes(k*(bg.η₀ - bg.η(x)))
            source = s_itp(x, k)
            s += sb * source * (xgrid[i+1] - xgrid[i])
        end
        return s
    end
# end

tEl_arr = zeros(length(camb_ℓ),length(ks_test))
for (i,ℓ) in enumerate(camb_ℓ)
    println(i)
    ℓð = sqrt((ℓ+2)*(ℓ+1)*ℓ*(ℓ-1)) #spin factor
    for (j,k) in enumerate(ks_test)
        bes = Bolt.bessel_interpolator(ℓ, k * bg_test.η₀)
        tEl_arr[i,j] = tEl(test_xi,k,all_test[8],bes,bg_test)
    end
end

#ℓ = 2
plot(camb_k, camb_E_ℓ_k[1,:],xscale=:log10,label="camb")
scatter!(ks_test, tEl_arr[1,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$E_{\ell=2}(k)$")
savefig("../../misc_plots/El_ℓ=2.png")

#ℓ = 40
plot(camb_k, camb_E_ℓ_k[26,:],xscale=:log10,label="camb")
scatter!(ks_test, tEl_arr[26,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$E_{\ell=40}(k)$")
savefig("../../misc_plots/El_ℓ=40.png")

#ℓ = 200
plot(camb_k, camb_E_ℓ_k[41,:],xscale=:log10,label="camb")
scatter!(ks_test, tEl_arr[41,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$E_{\ell=200}(k)$")
savefig("../../misc_plots/El_ℓ=200.png")

#ℓ = 2000
plot(camb_k, camb_E_ℓ_k[end-10,:],xscale=:log10,label="camb")
scatter!(ks_test, tEl_arr[end-10,:],xscale=:log10,label="bolt")
xlabel!(raw"$k$")
ylabel!(raw"$E_{\ell=2000}(k)$")
savefig("../../misc_plots/El_ℓ=2000.png")



# Structured runs
# the scheme is a list with [(nx_1,nk_1), (nx_1,nk_2),...,(nx_2,nk_1),...]
rrs = [single_experiment(nx,nk) 
        for nk in nkk 
         for nx in nxx]

# channge surface plot angle
plot(all_test[7](bg_test.x_grid,log10.(ks_test)))
surface(all_test[7](bg_test.x_grid,ks_test),camera=(130,20))
surface(all_test[8](bg_test.x_grid,ks_test),camera=(130,20))
surface(all_test[8](bg_test.x_grid[1:end-20],ks_test[20:end]),camera=(130,20))

        

plot(ℓs, camb_Cᵀᵀ)
plot!(ℓs, all_test[4],label="midpoint TT")
plot(ℓs, camb_Cᴱᴱ)
plot!(ℓs, all_test[6])
plot(ℓs, all_test[4]./camb_Cᵀᵀ)
plot!(ℓs, all_test[6]./camb_Cᴱᴱ)
plot!(ℓs, all_test[5]./camb_Cᵀᴱ)

nxx_2d = [nxx[i_k] for i_k in 1:length(nkk) for i_x in 1:length(nxx)]
nkk_2d = vcat([nkk for i_x in 1:length(nxx)]...)

# can take "TT","TE","EE" or any combination, but gets cluttered
# also penultimate argument can be used to look at raw or percent residual (but may need to adjust ylims)
save_path = "../../misc_plots/cl_conv_tests/"
p = plot_experiments(rrs,nxx_2d,nkk_2d,nxx,true,["TT"]) 
savefig(p,save_path*"rel_cltt_diffs.pdf")
savefig(p,save_path*"rel_cltt_diffs.png")
p = plot_experiments(rrs,nxx_2d,nkk_2d,nxx,true,["TE"]) 
savefig(p,save_path*"rel_clte_diffs.pdf")
savefig(p,save_path*"rel_clte_diffs.png")
p = plot_experiments(rrs,nxx_2d,nkk_2d,nxx,true,["EE"]) 
savefig(p,save_path*"rel_clee_diffs.pdf")
savefig(p,save_path*"rel_clee_diffs.png")

p = plot_experiments(rrs,nxx_2d,nkk_2d,nxx,false,["TT"]) 
savefig(p,save_path*"raw_cltt_diffs.pdf")
savefig(p,save_path*"raw_cltt_diffs.png")
p = plot_experiments(rrs,nxx_2d,nkk_2d,nxx,false,["TE"]) 
savefig(p,save_path*"raw_clte_diffs.pdf")
savefig(p,save_path*"raw_clte_diffs.png")
p = plot_experiments(rrs,nxx_2d,nkk_2d,nxx,false,["EE"]) 
savefig(p,save_path*"raw_clee_diffs.pdf")
savefig(p,save_path*"raw_clee_diffs.png")

# save the results
writedlm(save_path*"clttteee_diffs_v1.dat",rrs)



