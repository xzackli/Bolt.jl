using Bolt
using BenchmarkTools

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)
hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, 340bg.H₀)

xᵢ = log(1e-8)
u₀ = Bolt.initial_conditions(xᵢ, hierarchy)
udummy = deepcopy(u₀)
@btime Bolt.hierarchy!(udummy, u₀, hierarchy, xᵢ)

##
@btime  boltsolve(hierarchy);


##
k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
@time sf = source_grid(par, bg, ih, k_grid, BasicNewtonian());

## this is currently the expensive part
@time cltt(50:20:2000, par, bg, ih, sf);

##
