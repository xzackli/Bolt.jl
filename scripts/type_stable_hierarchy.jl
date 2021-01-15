using Bolt
# using PyPlot
# using Unitful, UnitfulAstro, NaturallyUnitful
# using Interpolations

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(SahaPeebles(), par, bg)
hierarchy = Hierarchy(340bg.H₀, par, bg, ih)

# using OrdinaryDiffEq
xᵢ = log(1e-8)
u₀ = Bolt.basic_newtonian_adiabatic_initial(hierarchy.par, hierarchy.bg, hierarchy.ih, xᵢ, 340bg.H₀)
udummy = deepcopy(u₀)
@btime Bolt.basic_newtonian_hierarchy!(udummy, u₀, hierarchy, xᵢ)
