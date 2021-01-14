using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful
using QuadGK


par = CosmoParams()
ℋ = x -> Bolt.ℋ(x, par)
ℋ′ = x -> Bolt.ℋ′(x, par)
# quadgk( ap -> 1.0 / (ap * Bolt.ℋ_a(ap, par)), 0.0, 1.0)


##
xx = -10:0.1:0.0
xmid = xx[1:end-1] .+ 0.05
clf()
plot(xmid, abs.(diff(ℋ.(xx)) ./ diff(xx)) )
plot(xmid, abs.(ℋ′.(xmid)), "-")
yscale("log")
gcf()
##
