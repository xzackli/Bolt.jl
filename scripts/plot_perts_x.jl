using Revise
using Bolt
using Plots
using Printf
using Interpolations

𝕡 = CosmoParams()
n_q=15
bg = Background(𝕡;x_grid=-20.0:0.1:0.0,nq=n_q)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
logqmin,logqmax = -6,-1
x_grid = collect(-10:0.1:0)

bg.ρ₀ℳ(-20.0)
bg.ρ₀ℳ(0.0)
bg.η(-20.0)
bg.ℋ(-20.0)

k/bg.ℋ(-9.02) #horizon entry
k/bg.ℋ(-5.4) /.7 #maybe what is happening is the fluid approximaiton turns on in class at thd eivergence?


Tγ = (15/ π^2 *bg.ρ_crit *𝕡.Ω_r)^(1/4)
νfac = (90 * 1.2020569 /(11 * π^4)) * (𝕡.Ω_r * 𝕡.h^2 / Tγ) *((𝕡.N_ν/3)^(3/4))
#^the factor that goes into nr approx to neutrino energy density, plus equal sharing ΔN_eff factor for single massive neutrino
Ω_νm = 𝕡.Σm_ν*νfac/𝕡.h^2
Ω_ν =  7*(2/3)*𝕡.N_ν/8 *(4/11)^(4/3) *𝕡.Ω_r
Ω_ν/Ω_νm/2 *(3/3.046)
abs.(class_pxs[6,:])[1]/abs.(ℳρ* 𝕡.h )[end]

ℓᵧ=100
ℓ_ν=100
ℓ_mν=50
reltol=1e-5 #cheaper  rtol
k =  1000bg.H₀*.3/.333 /10
kbolt = k/(bg.H₀*3e5/100)
xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))]

println("k = ", kbolt,
        " log10k = ", log10(kbolt), " h/Mpc")
pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5

results=zeros(pertlen,length(x_grid))
ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid))

hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
perturb = boltsolve(hierarchy; reltol=reltol)

for (i_x, x) in enumerate(x_grid)
    println(i_x)
    #the result of a mindless copy, obviously all the history is recorded at one k!
    #hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
    #perturb = boltsolve(hierarchy; reltol=reltol)
    u = perturb(x)  #z this can be optimized away, save timesteps at the grid!
    results[:,i_x] = u #z should use unpack somehow
    ℳρ[i_x],ℳσ[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
                            results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
                            bg,exp(x),𝕡)
    #normalization for plotting, divide by integral of just momentum measure
    ℳρ[i_x]=ℳρ[i_x] ./ (ρ_σ(ones(length(bg.quad_pts)),
                                   zeros(length(bg.quad_pts)),
                                   bg,exp(x),𝕡)[1] )


end
results

# #recompute this while testing
# ℳρtest,ℳσtest = zeros(length(x_grid)),zeros(length(x_grid))
# for (i_x, x) in enumerate(x_grid)
#     ℳρtest[i_x],ℳσtest[i_x] = ρ_σ(results[2(ℓᵧ+1)+(ℓ_ν+1)+1:2(ℓᵧ+1)+(ℓ_ν+1)+n_q,i_x],
#                             results[2(ℓᵧ+1)+(ℓ_ν+1)+2*n_q+1:2(ℓᵧ+1)+(ℓ_ν+1)+3*n_q,i_x],
#                             bg,exp(-20),𝕡)
# end

at(x) = bg.ρ₀ℳ(x)*((x2a(x)<1/100 ? 1 : 2))#/2 #*x2a(x)^4#*4
plot(x_grid, log10.( bg.ρ₀ℳ.(x_grid).*x2a.(x_grid).^4 ))
plot!(x_grid, log10.( bg.ρ₀ℳ.(x_grid).*x2a.(x_grid).^3 ),ls=:dash)
plot(x_grid, log10.(at.(x_grid)))
#plot(x_grid, log10.(8(bg.ρ₀ℳ.(x_grid) +bg.P₀ℳ.(x_grid))  ./ bg.ρ₀ℳ.(x_grid)),ls=:dash)
plot!(x_grid, log10.((bg.P₀ℳ.(x_grid)).*x2a.(x_grid).^4 *3),ls=:dot)
ρnr = bg.ρ₀ℳ.(x_grid)- 3bg.P₀ℳ.(x_grid)
ρr = 3bg.P₀ℳ.(x_grid)
plot!(x_grid, log10.( bg.ρ₀ℳ.(x_grid)))
plot!(x_grid, log10.( ρr ),ls=:dash)
plot!(x_grid, log10.( 2ρnr),ls=:dot)
plot!(x_grid, log10.( (ρnr+ρr)),ls=:dot)
plot!(x_grid, log10.( (ρr)),ls=:dot)
plot!(x_grid, log10.( (2ρnr+ρr)),ls=:dash)
plot!(x_grid, log10.( (ρnr+bg.ρ₀ℳ.(x_grid))),ls=:dot)
plot!(x_grid, log10.((2bg.ρ₀ℳ.(x_grid).-3bg.P₀ℳ.(x_grid))),ls=:dot)
xlims!(-6,-2)
ylims!(-10,-3)

#CLASS perturbations
#CLASS keys:
#['k (h/Mpc)', 'd_g', 'd_b', 'd_cdm', 'd_ur', 'd_ncdm[0]', 'd_tot',
#'phi', 'psi', 't_g', 't_b', 't_cdm', 't_ur', 't_ncdm[0]', 't_tot']
# ret = open("./test/data/class_px_kp3.dat","r") do datafile
ret = open("./test/data/class_px_kp03.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
retnf = open("../compare/class_px_kp03_nofluid.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end

#the second column is just a repeated k value, so remember it and delete col
kclass = ret[2][1]
class_pxs = transpose(reduce(hcat,ret[1:end .!= 2]))
class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))
class_pxs
println("kclass is ", kclass, " kbolt is ",kbolt)

#quick look at these - copying similar syntax from plot perts k
#skipping velocities this time just for simplicity
#matter δ
plot(class_pxs[1,:],log10.(abs.(class_pxs[4,:])),
     label=raw"$\delta_{c,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid,log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h),
      label=raw"$h \delta_{\rm{Bolt}}$",ls=:dash)

#baryon δ_b
plot!(class_pxs[1,:],log10.(abs.(class_pxs[3,:])),
    label=raw"$\delta_{b,\rm{CLASS}}$")
plot!(x_grid,log10.(abs.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]* 𝕡.h)),
    label=raw"$h \delta_{b,\rm{Bolt}}$",ls=:dash)

#throw in space metric Φ also
plot(class_pxs[1,:],log10.(abs.(class_pxs[8,:])),
    label=raw"$\Phi_{\rm{CLASS}}$")
plot!(x_grid, log10.(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h),
      label=raw"$h \Phi_{\rm{Bolt}}$",ls=:dash)
vline!([xhor],ls=:dot,c=:black,label="k/[2πℋ(x)h]=1")


xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt (NR) - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/nr_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

#massless neutrino monopole 𝒩0
plot(class_pxs[1,:],log10.(abs.(class_pxs[5,:])),
     label=raw"$\nu_{0,\rm{CLASS}}$",
     legend=:topleft)
plot!(x_grid, log10.(abs.(results[2(ℓᵧ+1)+1,:]* 𝕡.h*4)),
      label=raw"$4 h \nu_{0,\rm{Bolt}}$",ls=:dash)

#photon Θ0 monopole
plot(class_pxs[1,:],log10.(abs.(class_pxs[2,:])),
      label=raw"$\Theta_{0,\rm{CLASS}}$")
plot!(x_grid, log10.(abs.(results[1,:]* 𝕡.h*4)),
      label=raw"$4 h \Theta_{0,\rm{Bolt,8}}$",ls=:dash)

#massive neutrino monopole ℳ0
# plot(class_pxs[1,:],log10.(abs.(class_pxs[6,:])),
#     label=raw"$m\nu_{0,\rm{CLASS}}$",#)#,
#     legend=:topleft)

plot(class_pxsnf[1,:],log10.(abs.(class_pxsnf[6,:])),
    label=raw"$m\nu_{0,\rm{CLASS,nf}}$",
    legend=:topleft)#)#,
plot!(x_grid, log10.(abs.(ℳρ* 𝕡.h)),
    label=raw"$4 h m\nu_{0,\rm{3 Bolt}}$",ls=:dash)
    #ls=:dot)
vline!([xhor],ls=:dot,c=:black,label=raw"$k/(2\pi a H h)=1$")
vline!([log(1/100)],ls=:dot,c=:black)

# plot!(x_grid, log10.(abs.(ℳρtest* 𝕡.h *4)),
#     label=raw"$\mathrm{rel norm} 4 h m\nu_{0,\rm{Bolt}}$",ls=:dash)
#reverse the arrays because ow complains
itpnuclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[6,:][end:-1:1])
#drop the last bolt element because arrays are strangely aligned...
plot(x_grid[1:end-1], ((-ℳρ* 𝕡.h)[1:end-1]./itpnuclass.(x_grid[1:end-1]) ))
println(typeof(class_pxs[1,:]), ' ', typeof(class_pxs[6,:]))
hline!([1],ls=:dot,color=:black)

#check Phi, delta
itpphiclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[8,:][end:-1:1])
plot(x_grid[1:end-1], (results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]* 𝕡.h )[1:end-1]./itpphiclass.(x_grid[1:end-1]) )
hline!([1],ls=:dot,color=:black)

itpdelclass = LinearInterpolation(class_pxsnf[1,:][end:-1:1],class_pxsnf[4,:][end:-1:1])
plot(x_grid[1:end-1], -(results[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+2,:]* 𝕡.h )[1:end-1]./itpdelclass.(x_grid[1:end-1]) )
hline!([1],ls=:dot,color=:black)

xlabel!(raw"$x$")
ylabel!(raw"$\delta_{i}(x)$")
title!("Compare CLASS - Bolt (R) - k=$(@sprintf("%.3f", kclass))")
savefig("../compare/r_both_class_bolt_perts_x_k$(@sprintf("%.3f", kclass)).png")

minimum(class_pxs[1,:]),maximum(class_pxs[1,:])
minimum(x_grid),maximum(x_grid[1:end-1])


#utilities for mapping between comoving momenta and unit interval
from_ui(x,lqmi,lqma) = lqmi + (lqma- lqmi) / (1- (-1)) * (x- (-1))
xq2q(x,logqmin,logqmax) = 10.0 ^ from_ui(x,logqmin,logqmax)
Tν =  (𝕡.N_ν/3)^(1/4) *(4/11)^(1/3) * (15/ π^2 *bg.ρ_crit *𝕡.Ω_r)^(1/4)
logqmin,logqmax=log10(Tν/30),log10(Tν*30)
xq2q(10,logqmin,logqmax)
ftest = [dlnf0dlnq(xq2q(q,logqmin,logqmax),𝕡) for q in bg.quad_pts]
ai = exp(-20)
norm𝒩= 1 /( bg.ρ₀ℳ(log(ai))  *  ai^4)
normi = ρ_σ(ftest,
       zeros(length(bg.quad_pts)),
       bg,ai,𝕡)[1] * ai^-4 #* norm𝒩
onenormi = ρ_σ(ones(length(bg.quad_pts)),
       zeros(length(bg.quad_pts)),
       bg,ai,𝕡)[1]* ai^-4 #* norm𝒩

normi/onenormi

normf= ρ_σ(ftest,
       zeros(length(bg.quad_pts)),
       bg,1,𝕡)[1]
onenormf = ρ_σ(ones(length(bg.quad_pts)),
       zeros(length(bg.quad_pts)),
       bg,1,𝕡)[1]

normf/onenormf


#consistency check on background at initial and final times
Ω_ν*ai^(-4) *bg.ρ_crit
bgi = bg.ρ₀ℳ(log(ai)) #extra factor for 2 massless neutrinos

Ω_νm*bg.ρ_crit
bgf = bg.ρ₀ℳ(0)

normi / bgi / 4
normf / bgf /4

normi / ((2bg.ρ₀ℳ(log(ai))-3bg.P₀ℳ(log(ai)))) / 4
normf / ((2bg.ρ₀ℳ(0)-3bg.P₀ℳ(0))) /4

norm𝒩 =  ρ_σ(map(dlnf0dlnq,(bg.quad_pts,𝕡),zeros(length(bg.quad_pts)),bg,ai,𝕡)[1] * a^-4
