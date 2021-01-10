using Bolt
using PyPlot
using NLsolve
using Unitful, UnitfulAstro, NaturallyUnitful

##
par = Cosmo()
z_grid = 1587.4:10.0:1800
a_grid = 1.0 ./ (z_grid .+ 1)
clf()
plot( z_grid, Bolt.saha_Xₑ(a_grid, par), "-", label="Saha")
ylabel(raw"$X_e$")
legend()
gcf()

##

Bolt.peebles_rhs(0.99, 1/1587, par)


##
Bolt.peebles_Xₑ(par, 0.99, 1/1587, 0.2)

##
using DifferentialEquations
f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
