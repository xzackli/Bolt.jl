using Bolt
using PyPlot
using Unitful, UnitfulAstro, NaturallyUnitful

import Bolt: z2a

par = Cosmo()

# integrate saha to some transition redshift, we choose 1587.4
z_transition = 1587.4
a_transition = z2a(z_transition)
saha_z_grid = 1800:-10:z_transition
peebles_z_grid = z_transition:-10:100

early_time_Xₑ = Bolt.saha_Xₑ(par)
late_time_Xₑ = Bolt.peebles_Xₑ(
    par, early_time_Xₑ(a_transition), a_transition, 0.01)

clf()
plot(saha_z_grid, [early_time_Xₑ(z2a(z)) for z in saha_z_grid], "-", label="Saha")
plot(peebles_z_grid, [late_time_Xₑ(z2a(z)) for z in peebles_z_grid], "-", label="Peebles")
ylabel(raw"$X_e$")
xlim(1800,0)
ylim(2e-4, 2)
yscale("log")
xlabel(raw"$z$")
legend()
gcf()
