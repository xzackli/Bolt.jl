using Bolt
using BenchmarkTools
using OrdinaryDiffEq

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)

k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
k_i = 1
hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k_grid[k_i])

# xᵢ = log(1e-7)
xᵢ = (hierarchy.bg.x_grid)[10]
u₀ = Bolt.initial_conditions(xᵢ, hierarchy)
prob = ODEProblem{true}(Bolt.hierarchy!, u₀, (xᵢ, 0.0), hierarchy)
#@time
sol = solve(prob, KenCarp3(), reltol=1e-4)
