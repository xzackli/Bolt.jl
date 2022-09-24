using Bolt
using ForwardDiff
using Plots
using BenchmarkTools
using Plots.PlotMeasures


# General derivates of  all other parameters
function clp(𝕡::DT,ells) where DT
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=15)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    return cltt(ells, 𝕡, bg, ih, sf)
end


pt = CosmoParams() #dummy struct
# do a loop over everyone to see who is slow...this is bad for performance but who cares we are just checking times inside...
fds = zeros(length(fieldnames(typeof(pt))),length(ells))
fws = zeros(length(fieldnames(typeof(pt))),length(ells))

for (i,name) in enumerate(fieldnames(typeof(pt)))
    println("i = ", i, " $name = ",getfield(pt,name))
    xv = getfield(pt,name)#pt.name
    function clx(x::DT,name,ells) where DT
        𝕡n = CosmoParams{DT}(;Dict((fieldnames(typeof(pt))[i])=>x)...)
        return clp(𝕡n,ells)
    end
    f(x) = clx(x,name,ells)
    fwddiff_∂cl = ForwardDiff.derivative(f,xv)
    println(fwddiff_∂cl[1])
    Δ = 1e-2*xv 
    finitediff_∂cl = (f(xv + Δ) .- f(xv - Δ)) ./ 2Δ
    fws[i,:] = fwddiff_∂cl
    fds[i,:] = finitediff_∂cl
end

#multiply the mnu by its value to get something of OOM as the others...
#just an ugliness of the Mpc units for the eV neutrino mass...
fws[end,:] .*= 𝕡.Σm_ν
fds[end,:] .*= 𝕡.Σm_ν

# save the results
writedlm("ClTT_AD.dat",fws)
writedlm("ClTT_FD1e-2.dat",fds)

#make some plots
pnames = string.( fieldnames(typeof(pt)) )
for i in 1:length(fieldnames(typeof(pt)))
    plot(title=raw"$ p = $"*pnames[i])
    plot!(ells, fws[i,:] .* ells.^2,label="AD")
    plot!(ells, fds[i,:] .* ells.^2,ls=:dash,label="FD, "*raw"$\frac{\Delta p}{p} = 10^{-2}$")
    xlabel!(raw"$\ell$")
    ylabel!(raw"$\ell^{2} \frac{\partial C_{\ell}^{TT}}{\partial p }$",left_margin=5mm)
    savefig("AD_vs_FD1e-2_all_params_"*pnames[i]*".pdf")
end



#-----FIXME move the below, which just tests scalar polzn spectra, and has nothing to do with derivs

𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=15)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
ih = IonizationHistory(𝕣, 𝕡, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
sf_t = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
sf_e = source_grid_P(𝕡, bg, ih, k_grid, BasicNewtonian())

ℓs = 10:20:1200
Cℓtt = cltt(ℓs, 𝕡, bg, ih, sf_t)
Cℓte = clte(ℓs, 𝕡, bg, ih, sf_t,sf_e)
Cℓee = clee(ℓs, 𝕡, bg, ih, sf_e)

ℓfac = ℓs.*(ℓs.+1)
plot(ℓs, @. ( ℓfac * Cℓtt))
xlabel!(raw"$\ell$")
ylabel!(raw"$\ell(\ell+1)C^{TT}_{\ell}$")
# savefig("./test/cltt.png")
plot(ℓs, @. (ℓfac *Cℓte))
xlabel!(raw"$\ell$")
ylabel!(raw"$\ell(\ell+1)C^{TE}_{\ell}$")
# savefig("./test/clte.png")
plot(ℓs, @. (ℓfac * Cℓee))
xlabel!(raw"$\ell$")
ylabel!(raw"$\ell(\ell+1)C^{EE}_{\ell}$")
# savefig("./test/clee.png")


#--- old code
# Cₗ as a function of baryon density
#Some noise at lowest few ells...9/17/21

function clb(Ω_b::DT, ells) where DT
    𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=15)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    sf = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    return cltt(ells, 𝕡, bg, ih, sf)
end

ells = 10:20:1200
f(Ω_b) = clb(Ω_b, ells)
@time cl = f(0.046)
@time ∂cl = ForwardDiff.derivative(f, 0.046)  # you can just ForwardDiff the whole thing

@profview f(0.046)
#58.130749 seconds (74.01 M allocations: 7.458 GiB, 3.10% gc time, 2.92% compilation time)
@profview ForwardDiff.derivative(f, 0.046)
#219.854089 seconds (139.67 M allocations: 13.256 GiB, 1.45% gc time, 8.75% compilation time)


ells = 10:20:1200
fh(h) = clh(h, ells)
@time clhh = fh(0.7)
#^ 59.244405 seconds (72.37 M allocations: 7.376 GiB, 3.04% gc time, 0.28% compilation time)
@time ∂clhh = ForwardDiff.derivative(fh, 0.7)  # you can just ForwardDiff the whole thing
#^ 217.447492 seconds (109.88 M allocations: 12.040 GiB, 1.53% gc time, 9.47% compilation time)

@profview fh(0.7)
@profview ForwardDiff.derivative(fh, 0.7)
