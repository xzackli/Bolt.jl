#jms - Marius's notebook
# setup
using Revise
using Bolt
# using Plots
# using Interpolations
using DelimitedFiles
using FFTW
using PyCall
using LinearAlgebra
using Roots

#bolt setup
𝕡 = CosmoParams()
bg = Background(𝕡)#; x_grid=-25.0:0.01:0.0)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
# x_grid = -25:0.01:0.0

#time and space
nk = 1000
nx = 200
z_span = (10000000,10000)

#runs from -16.1 -> -9ish
x_grid = map(range(bg.η.(z2x.(z_span))..., length=100)) do η
	find_zero(x -> bg.η.(x) - η, z2x.(z_span))
end
#big k grid
k_grid = range(0, 20000bg.H₀, length=(nk+1))[2:end]
k_grid_full = [-reverse(k_grid); 0; k_grid]
#what does this @ do? oh cool in does the dot at every point - convenient!
Φx₀ = @. exp(-(k_grid_full/$last(k_grid))^2/(2*0.01^2))
Φk₀ = rfft(Φx₀)[2:end]

solns = map(enumerate(k_grid), Φk₀) do (i,k), Φ₀
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, 8)
    perturb = boltsolve(hierarchy; reltol=1e-8, Φ₀=abs(Φ₀))
	(perturb,hierarchy)
end;

#some thing to do the FT and unpack all at once?
#change it to account for ℳ
using Bolt: unpack
grid = broadcast(solns, Φk₀, x_grid') do (perturb,hierarchy), Φ₀, x
	Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = unpack(perturb(x), hierarchy)
	exp(im*angle(Φ₀)) .* ((Θ[0] - 0), (𝒩[0] - 0), δ - 0, δ_b - 0)
end;

plot(grid[:,1])
plot(grid[:,])



plot(grid_real_θ[:,80] ./(5*minimum(grid_real_θ[:,80])))
plot!(grid_real_δ_b[:,80] ./(5*minimum(grid_real_δ_b[:,80])))
ylims!(-.1,.3)
xlims!(1000,1200)
#look
zz = [ones(100); zeros(nk-100+1)]
y = [zero(x_grid)'; getindex.(grid,1)][:,end] .* zz
using Plots
plot(-irfft(y,2nk+1))

#pull out components we want to look at
grid_real_θ   = irfft([zero(x_grid)'; getindex.(grid,1)] .* zz, 2nk+1, 1);
grid_real_ν   = irfft([zero(x_grid)'; getindex.(grid,2)],       2nk+1, 1);
grid_real_δ   = irfft([zero(x_grid)'; getindex.(grid,3)],       2nk+1, 1);
grid_real_δ_b = irfft([zero(x_grid)'; getindex.(grid,4)] .* zz, 2nk+1, 1);

#save for animating later
writedlm("../compare/dcosmoviz_theta.dat",grid_real_θ)
writedlm("../compare/dcosmoviz_nu.dat",grid_real_ν)
writedlm("../compare/dcosmoviz_delta.dat",grid_real_δ)
writedlm("../compare/dcosmoviz_deltab.dat",grid_real_δ_b)
writedlm("../compare/dcosmoviz_xgrid.dat",x_grid)
