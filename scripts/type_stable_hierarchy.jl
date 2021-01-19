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
