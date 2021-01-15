using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using Interpolations

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(SahaPeebles(), par, bg)
hierarchy = Hierarchy(340bg.H₀, par, bg, ih)

# @time sol = solve(prob, Rodas5(), reltol=1e-10)
S = solve(hierarchy, BasicNewtonian(); reltol=1e-10)
