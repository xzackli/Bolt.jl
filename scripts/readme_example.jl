#readme example without diff, just want the Cltt
using Revise
using Bolt
using Plots
using Interpolations
using DelimitedFiles

#modified version of readme example - looks good

function source_grid_t(par::AbstractCosmoParams{T}, bg, ih, k_grid; ℓᵧ=10, reltol=1e-11) where T
    x_grid = bg.x_grid
    grid = zeros(T, length(x_grid), length(k_grid))
    for (i_k, k) in enumerate(k_grid)
        println("i_k = ", i_k)
        hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k, ℓᵧ)
        perturb = boltsolve(hierarchy; reltol=reltol)
        for (i_x, x) in enumerate(x_grid)
            u = perturb(x)  # this can be optimized away, save timesteps at the grid!
            du = similar(u)
            Bolt.hierarchy!(du, u, hierarchy, x)
            grid[i_x,i_k] = Bolt.source_function(du, u, hierarchy, x)
        end
    end
    # return grid
    itp = LinearInterpolation((x_grid, k_grid), grid, extrapolation_bc = Line())
    return itp
end

𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-25.0:0.1:0.0)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
sf = source_grid_t(𝕡, bg, ih, k_grid)
ℓs = 0:1:5200
cells = [cltt(ℓ, 𝕡, bg, ih, sf) for ℓ in ℓs]

# plot and save
plot(log10.(ℓs),ℓs.^2 .* cells)
xlims!(0,3000)
writedlm("../compare/5200_bolt_cmb_cell.txt", (ℓs,cells))



#playing around below

bg.η(bg.x_grid[end])
x2η = spline( [bg.η(x) for x in bg.x_grid], bg.x_grid)
function η2a(η,dx=0.01)
  x_fine = minimum(bg.x_grid):dx:maximum(bg.x_grid)
  return exp(x_fine[argmin(abs.(x2η.(x_fine).-η))])
end
a2η(a) = bg.η(log(a))


#checking the sound horizon...
function rs(η)
  xq,wq = bg.quad_pts,bg.quad_wts
  logηmin,logηmax=log10(x2η(minimum(bg.x_grid)))-10,log10(η)
  #η = ∫ (aℋ)^-1 da/dx
  Ics(y) = ih.csb²( log(η2a(xq2q(y,logηmin,logηmax))) )^(1/2)  / dxdq(xq2q(y,logηmin,logηmax),logηmin,logηmax)
  return sum(Ics.(xq).*wq)
end

plot!(bg.x_grid,log10(.75*π) .+ log10.(bg.η(maximum(bg.x_grid))) .- log10.(rs.(x2η.(bg.x_grid))))

rs(x2η(log(1/(1+1100))))

plot(bg.x_grid,log10.((ih.csb²(bg.x_grid).^(1/2).*bg.ℋ.(bg.x_grid)).^(-1/2)))


plot(bg.x_grid,log10.(exp.(bg.x_grid).*(ih.csb²(bg.x_grid).^(1/2))))
vline!([log(1/(1+1100))])

#rS in eV
hundredkmsMpc_to_eV = 2.13e-33
kmsMpc_to_eV = hundredkmsMpc_to_eV/100
rs_try = rs(x2η(log(1/(1+1100)))) * kmsMpc_to_eV *3e5^2 /100  / 𝕡.h^2
rsMpch = rs(x2η(log(1/(1+1100))))*(bg.H₀*3e5/100)

#to convert from eV to our units
k_grid_hMpc = k_grid/(bg.H₀*3e5/100) #denom factor is equal to h * ev2kmsmpc100 * 3e5 /100 = h* ev2kmsMpc * 3e5
ℓpapprox=220 #not sure what exactly this should be
kpapprox = ℓpapprox/bg.η₀/.75 #dodelson pg 247
kpapprox_hMpc = kpapprox/(bg.H₀*3e5/100)
rsp1approx_hMpc = π/kpapprox_hMpc
rsp1approx_Mpc = rsp1approx_hMpc / .7
#this is at least the right ballpark ~140-150 Mpc - I don't expect to get it exactly bc using diff cosmo

#going backwards to check
ℓpapprox2 = 0.75*π*bg.η₀*bg.H₀*3e5/(100*110)

ℓpapprox3 = 0.75*π*bg.η₀*bg.H₀*3e5/(100*rs_try)

bg.η₀*(bg.H₀*3e5/100)  / .7#these are the class units of τ in Mpc for conf time
