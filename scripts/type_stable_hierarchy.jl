using Bolt
using BenchmarkTools

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)

k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
k_i = 10
hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k_grid[k_i])

xᵢ = log(1e-8)
u₀ = Bolt.initial_conditions(xᵢ, hierarchy)
udummy = deepcopy(u₀)
@btime Bolt.hierarchy!(udummy, u₀, hierarchy, xᵢ)

##
@time perturb = boltsolve(hierarchy);

##
@time sf = source_grid(par, bg, ih, k_grid, BasicNewtonian());

##
