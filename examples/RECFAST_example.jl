using Bolt
𝕡 = CosmoParams()
bg = Background(𝕡)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)  #  𝕣 = Bolt.Peebles()
ih = IonizationHistory(𝕣, 𝕡, bg)

dz = 10000/1000
z = (10000 - dz):(-dz):0.0

using Plots, LaTeXStrings
plot(z, ih.Tmat.(z2x.(z)), ls=:solid, label=L"T_{\mathrm{mat}}", xlim=(10,10000),
    xscale=:log10, xlabel="redshift", ylabel="temperature [K]", yscale=:log10, legend=:bottomright)
plot!(z, ih.Trad.(z2x.(z)), ls=:dash, label=L"T_{\mathrm{rad}}")
