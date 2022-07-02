
using Adapt, BenchmarkTools, Bolt, CUDA, ComponentArrays, Interpolations, 
    OffsetArrays, OrdinaryDiffEq, PyPlot, Setfield, StaticArrays

CUDA.allowscalar(false)

##

T = Float64
storage = CUDA.CuArray{T} # for GPU
# storage = nothing # for CPU

##

# bg/ion setup
𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=5)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k = 500bg.H₀
reltol=1e-5
hierarchy = adapt(storage, Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, 10, 10, 8, 5));

##

xᵢ = first(hierarchy.bg.x_grid)
u₀ = CUDA.@allowscalar adapt(storage, Bolt.initial_conditions(xᵢ, hierarchy));
du = zero(u₀);

##

cudacall_hierarchy!(du, u, h, x) = @cuda threads=1 blocks=1 Bolt.hierarchy!(du, u, h, x)

# works
@btime CUDA.@sync cudacall_hierarchy!(du, u₀, hierarchy, xᵢ)

##

prob = ODEProblem{true}(cudacall_hierarchy!, u₀, (T(xᵢ), T(0)), hierarchy);
# broken, probably need to define a few more ::ComponentArray{<:CuArray} stuff to make it work
sol = @time solve(prob, KenCarp4(), reltol=1e-4, saveat=hierarchy.bg.x_grid, dense=false)
