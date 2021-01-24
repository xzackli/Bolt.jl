using Bolt
using BenchmarkTools

par = CosmoParams()
bg = Background(par)
ih = IonizationHistory(Peebles(), par, bg)

k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
k_i = 10
hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k_grid[k_i])

xᵢ = log(1e-8)
u₀ = Bolt.initial_conditions(xᵢ, hierarchy)
udummy = deepcopy(u₀)
@btime Bolt.hierarchy!(udummy, u₀, hierarchy, xᵢ)

##
@time perturb = boltsolve(hierarchy);

##
@time sf = source_grid(par, bg, ih, k_grid, BasicNewtonian());

##

using Hankel
R = bg.η₀ * bg.H₀
N = 4096
q = QDHT{1,1}(R, N)
using DataInterpolations

x_from_ηH₀ = CubicSpline(bg.x_grid, bg.η.(bg.x_grid) .* bg.H₀)

function f(ηH₀)
    x = x_from_ηH₀(ηH₀)
    u = perturb(x)  # this can be optimized away, save timesteps at the grid!
    du = similar(u)
    Bolt.hierarchy!(du, u, hierarchy, x)
    return Bolt.source_function(du, u, hierarchy, x)
end
fr = f.(q.r)
# fk = q * fr # 0th-order QDHT => fk should have only the first entry non-zero
# fk .*= 1.1

##
clf()
# figure()

# xx =  bg.η.(bg.x_grid[1200:end-1]) .* bg.H₀

xx =  bg.η.(bg.x_grid[1200:end-1]) .* bg.H₀
plot(xx, abs.(f.(xx)) ./ xx,
    label=raw"$k/H_0=$" * "$(k_grid[k_i]/bg.H₀)", "-")

xx = q.r
# plot(xx, abs.(f.(xx)) ./ xx, ".",
#     label=raw"$k/H_0=$" * "$(k_grid[k_i]/bg.H₀)")

xlabel(raw"$η H_0$")
ylabel("Source Function")
xscale("log")
# xlim(0.05, 0.1)
legend()
# yscale("log")
# ylim(1e-5,800)
ylim(1e-5,20)
gcf()


##
bg.η₀ * k_grid[end]

##
clf()
figure()

# xx =  bg.η(bg.x_grid[1200:end]) .* bg.H₀
xx
plot(xx, abs.(sf[1200:end,k_i]) ./ xx, "-",
    label=raw"$k/H_0=$" * "$(k_grid[k_i]/bg.H₀)")
xlabel(raw"$η H_0$")
ylabel("Source Function")
# xscale("log")
xlim(0.05, 0.2)
legend()
# yscale("log")
# ylim(1e-5,20)
gcf()

##
## this is currently the expensive part
# function test(par, bg, ih, sf)
#     @time cltt(2000, par, bg, ih, sf)
# end

# test(par, bg, ih, sf)
