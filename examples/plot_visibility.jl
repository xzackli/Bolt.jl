using Bolt
using Plots, LaTeXStrings

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)
x_grid = bg.x_grid

l = @layout [a ; b]
p1 = plot(x_grid, ih.τ.(x_grid), label=L"\tau", yscale=:log10)
plot!(x_grid, abs.(ih.τ′.(x_grid)), ls=:dash, label=L"|\tau^\prime|")
p2 = plot(x_grid, ih.g̃.(x_grid), label=raw"$\tilde{g}$", xlim=(-8.0, -6.0), xlabel="x")
plot!(x_grid, ih.g̃′.(x_grid) ./ 10, ls=:dash, label=raw"$\tilde{g}\prime/10$")
plot!(x_grid, ih.g̃′′.(x_grid) ./ 300, ls=:dash, label=raw"$\tilde{g}\prime/300$")
plot(p1, p2, layout = l, size=(350,500))
