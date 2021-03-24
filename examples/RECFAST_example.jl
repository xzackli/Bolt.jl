using Bolt
𝕡 = CosmoParams()
bg = Background(𝕡)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)  #  𝕣 = Bolt.Peebles()
ih = IonizationHistory(𝕣, 𝕡, bg)

Nz = 1000
dz = 10000/float(Nz)
z = (10000 - dz):(-dz):0.0
##

using PyPlot
clf()
plot(z, ih.Tmat.(z2x.(z)), "-", label=raw"$T_{\mathrm{mat}}$")
plot(z, ih.Trad.(z2x.(z)), "--", label=raw"$T_{\mathrm{rad}}$")

yscale("log")
xscale("log")
legend()
ylim(1, 2e4)
xlim(10, 10000)
xlabel("redshift")
ylabel("temperature [K]")
gcf()
