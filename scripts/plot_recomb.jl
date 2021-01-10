using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful

par = Cosmo()

# integrate saha to some transition redshift, we choose 1587.4
z_transition = 1587.4
saha_z_grid = 1800:-10:z_transition
saha_a_grid = 1.0 ./ (saha_z_grid .+ 1)
saha_Xₑ_grid = Bolt.saha_Xₑ(saha_a_grid, par)

peebles_a_grid, peebles_Xₑ_grid = Bolt.peebles_Xₑ(
    par, saha_Xₑ_grid[end], saha_a_grid[end], 0.01)

clf()
plot(saha_z_grid, saha_Xₑ_grid, "-", label="Saha")
plot(1 ./ peebles_a_grid .- 1, peebles_Xₑ_grid, "-", label="Peebles")
ylabel(raw"$X_e$")
xlim(1800,0)
ylim(2e-4, 2)
yscale("log")
xlabel(raw"$z$")
legend()
gcf()
