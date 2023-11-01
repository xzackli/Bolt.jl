# using Revise
using Bolt
using ForwardDiff
using Plots
using BenchmarkTools
using Printf
using Interpolations,DataInterpolations
using Plots.PlotMeasures

#input ingredients
𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)  #  𝕣 = Bolt.Peebles()
ih = IonizationHistory(𝕣, 𝕡, bg)
logqmin,logqmax = -6,-1
logq_pts = logqmin:(logqmax-logqmin)/(n_q-1):logqmax

kmin,kmax= 0.1bg.H₀*100,5000bg.H₀
k_grid = log10_k(kmin,kmax,33)

ℓᵧ=25 #cutoff
ℓ_ν=25#10
ℓ_mν=20
reltol=1e-8 #cheaper  rtol
x=0
a=exp(x)
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
results=zeros(pertlen,length(k_grid))
for (i_k, k) in enumerate(k_grid)
    println(i_k)
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
    perturb = boltsolve(hierarchy; reltol=reltol)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_k] = u #z should use unpack somehow
    # results[:,i_k] .= boltsolve_rsa(hierarchy;reltol=reltol)[:,end]
end
results

#Integrate the q moments for ℳ0 and ℳ1 for gauge transformatoin/plotting
ℳρ,ℳθ = zeros(length(k_grid)),zeros(length(k_grid))
for (i_k, k) in enumerate(k_grid)
    ℳρ[i_k],_ = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_k],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_k],
                            bg,a,𝕡)./ bg.ρ₀ℳ(x)
    #Below assumes negligible neutrino pressure for the normalization (fine at z=0)
    ℳθ[i_k] = k*θ(results[2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+2n_q,i_k],
                     bg,a,𝕡)./ bg.ρ₀ℳ(x)
    #Also using the fact that a=1 at z=0
end

k_grid_hMpc = k_grid/(bg.H₀*299792.458/100)#3e5/100)
#^convert our units of eV/(eV 100 h /c) to km/s/Mpc/h

#put together the matter transfer function (Newtonian gauge)
#Newtonian perturbations
δcN,δbN = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]
vcN,vbN = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+3,:],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5,:]
ℳρN,ℳθN = ℳρ,ℳθ
vmνN = -ℳθN./ k_grid

#omegas to get weighted sum for total matter in background
Tγ = (15/ π^2 *bg.ρ_crit *𝕡.Ω_r)^(1/4)
ζ = 1.2020569
νfac = (90 * ζ /(11 * π^4)) * (𝕡.Ω_r * 𝕡.h^2 / Tγ) *((𝕡.N_ν/3)^(3/4))
#^the factor that goes into nr approx to neutrino energy density, plus equal sharing ΔN_eff factor for single massive neutrino
Ω_ν = 𝕡.Σm_ν*νfac/𝕡.h^2
Ωm = 𝕡.Ω_m+𝕡.Ω_b+Ω_ν

#construct gauge-invariant versions of density perturbations
δc = δcN - 3bg.ℋ(x)*vcN./k_grid
δb = δbN - 3bg.ℋ(x)*vbN./k_grid
#assume neutrinos fully non-relativistic and can be described by fluid (ok at z=0)
δmν = ℳρN - 3bg.ℋ(x)*vmνN./k_grid
#FIXME: what does this mean as the fluid approximation gets worse?

#checking what it looks like when doing gauge change
#b
plot(log10.(k_grid_hMpc),log10.(δbN))
plot!(log10.(k_grid_hMpc),log10.(3bg.ℋ(x)*vbN./k_grid),ls=:dot)
plot!(log10.(k_grid_hMpc),log10.(abs.(δb)),ls=:dot)


#c
plot!(log10.(k_grid_hMpc),log10.(δcN))
plot!(log10.(k_grid_hMpc),log10.(3bg.ℋ(x)*vcN./k_grid),ls=:dot)
plot!(log10.(k_grid_hMpc),log10.(abs.(δc)),ls=:dot)
#ν
plot!(log10.(k_grid_hMpc),log10.(ℳρN))
plot!(log10.(k_grid_hMpc),log10.(3bg.ℋ(x)*vmνN./k_grid),ls=:dot)
plot!(log10.(k_grid_hMpc),log10.(abs.(δmν)),ls=:dot)

#check that gauge transformation results in 1 at sub-horizon scales
plot(log10.(k_grid_hMpc),k_grid_hMpc.*δbN ./ δb,label="b")
ylims!(0.0008,.002)
plot!(log10.(k_grid_hMpc),k_grid_hMpc.*δcN ./ δc,label="c")
plot!(log10.(k_grid_hMpc),k_grid_hMpc.*ℳρN ./ δmν,label="n")
hline!([1],c="black",ls=:dot)

#put together gauge-invariant matter
δm = (𝕡.Ω_m*δc + 𝕡.Ω_b*δb + Ω_ν*δmν) ./ Ωm
As=1e-10*exp(3.043)
Pprim = As*(k_grid_hMpc./0.05).^(𝕡.n-1)
PL_un = (2π^2 ./ k_grid_hMpc.^3).*(δm).^2 .*Pprim

#Load and plot the CLASS matter power at z=0
ret = open("/Users/jsull/Documents/berkeley/bolt/compare/zack_N_class_pk_dense_x0_nofluid_nonu.dat","r") do datafile
# ret = open("../compare/class_pk_x0_nofluid.dat","r") do datafile
# ret = open("../compare/class_pk_xm5_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
class_pk = reduce(hcat,ret)

#plot power
itpclass = CubicSpline(class_pk[2,:], class_pk[1,:])
p1 = plot(log10.(k_grid_hMpc[minkcut:maxkcut]),log10.(PL_un[minkcut:maxkcut]),label="Bolt")#,marker=:dot)
# plot!(log10.(class_pk[1,:]),log10.(class_pk[2,:]),label="CLASS",ls=:dash)
plot!(log10.(k_grid_hMpc[minkcut:maxkcut]),log10.(itpclass.(k_grid_hMpc[minkcut:maxkcut])),label="CLASS",ls=:dash)
ylabel!(raw"$\log ~P_{L}(k)$")

#interpolate to get ratio
minkcut=1
maxkcut=length(k_grid_hMpc)-3
p2=plot(log10.(k_grid_hMpc[minkcut:maxkcut]), (PL_un[minkcut:maxkcut])./itpclass.(k_grid_hMpc[minkcut:maxkcut]),
        legend=false,left_margin=4mm)#,marker=:circle )
hline!([1],ls=:dash,c=:black)
hline!([1.0001],ls=:dot,c=:black)
hline!([0.9999],ls=:dot,c=:black)
ylims!(.999,1.001)
ylabel!(raw"$\frac{P_{L,\rm{Bolt}}}{P_{L,\rm{CLASS}}}(k)$")
xlabel!(raw"$\log ~k \ [h/Mpc]$")
# xlims!(log10(minimum(class_pk[1,:])),log10(maximum(class_pk[1,:])))
# (2*PL_un[minkcut:end])./itpclass.(k_grid_hMpc[minkcut:end])

l = @layout [a  ; b]
plot(p1, p2, layout = l)
title!("Plin CLASS - Bolt - z=$(@sprintf("%.0f", exp(-x)-1))")
savefig("/Users/jsull/Documents/berkeley/bolt/compare/ellmax25_plin_both_class_bolt_perts_k_z$(@sprintf("%.0f", exp(-x)-1)).png")

#----
#This is old
#WIP try ForwardDiff copying Cl_TT code
#Derivative of Ωm for plot

#need to streamline this, should have As as a Bolt parameter (and eventually σ8)
#generalize the perturbations to be computed from hierarchy
function plin(k,results)
    #copy code abvoe
    ℳρ,_ = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q],
                            bg,a,𝕡)./ bg.ρ₀ℳ(x)
    #Below assumes negligible neutrino pressure for the normalization (fine at z=0)
    ℳθ = k*θ(results[2(ℓᵧ+1)+(ℓ_ν+1)+n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+2n_q],
                     bg,a,𝕡)./ bg.ρ₀ℳ(x)
    #Also using the fact that a=1 at z=0
    δcN,δbN = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h
    vcN,vbN = results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+3,:],results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5,:]* 𝕡.h
    ℳρN,ℳθN = ℳρ,ℳθ
    vmνN = -ℳθN./ k
    #omegas to get weighted sum for total matter in background
    Tγ = (15/ π^2 *bg.ρ_crit *𝕡.Ω_r)^(1/4)
    ζ = 1.2020569
    νfac = (90 * ζ /(11 * π^4)) * (𝕡.Ω_r * 𝕡.h^2 / Tγ) *((𝕡.N_ν/3)^(3/4))
    #^the factor that goes into nr approx to neutrino energy density, plus equal sharing ΔN_eff factor for single massive neutrino
    Ω_ν = 𝕡.Σm_ν*νfac/𝕡.h^2
    Ωm = 𝕡.Ω_m+𝕡.Ω_b+Ω_ν

    #construct gauge-invariant versions of density perturbations
    δc = δcN - 3bg.ℋ(x)*vcN./k
    δb = δbN - 3bg.ℋ(x)*vbN./k
    #assume neutrinos fully non-relativistic and can be described by fluid (ok at z=0)
    δmν = ℳρN - 3bg.ℋ(x)*vmνN./k
    println(δb,δmν)
    δm = (𝕡.Ω_m*δc .+ 𝕡.Ω_b*δb .+ Ω_ν*δmν) ./ Ωm
    As=1e-10*exp(3.043)
    k_hMpc=k/(bg.H₀*3e5/100)
    Pprim = As*(k_hMpc./0.05).^(𝕡.n-1)
    PL= (2π^2 ./ k_hMpc.^3).*(δm*𝕡.h).^2 .*Pprim
    return PL
end

#PL as a function of Ωmf
#one k at a time to start
function PL(Ω_c::DT, k) where DT
   𝕡 = CosmoParams{DT}(Ω_m=Ω_c)
   println(𝕡)
   bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=n_q)
   𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
   ih = IonizationHistory(𝕣, 𝕡, bg)
   #Why does recfast not work?
   # ih = IonizationHistory(Peebles(), 𝕡, bg)
   #Peebles doesn't work for same reason as before, raw Ha(a) function,
   #Can't fix this easily either because function that takes it for Peebles
   #doesn't know about bg...
   hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
   perturb = boltsolve(hierarchy; reltol=reltol)
   u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
   results=zeros(pertlen)
   results = u
   # println(results)
   return plin(k, results)
end

f(Ω_c) = PL(Ω_c, k_grid[1])#
#this is gonna take forever
@time pl = f(0.224)
@time ∂pl = ForwardDiff.derivative(f, 0.224)
