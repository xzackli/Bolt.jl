using Quadrature, Cubature, Zygote, ForwardDiff
using Test
using Bolt


xgrid = 0.0:0.001:π

function approx_cos(x, p, xgrid)
    y = [cos(p * x_) for x_ in xgrid]
    f = Bolt.spline(y, xgrid)
    f(x)
end

##
g = x -> approx_cos(x, 2.0, xgrid)

g'(π/4)

# clf()
# plot(xgrid, -sin.(xgrid), "-")
# plot(xgrid, g'.(xgrid), "--")
# gcf()

##




### One Dimensional

# par = CosmoParams()

# lb = 1.0
# ub = 3.0
# p = 2.0
# # prob = QuadratureProblem(f,lb,ub,p)
# # sol = solve(prob,QuadGKJL(),reltol=1e-3,abstol=1e-3)

# function testf(x)
#     res = Zygote.Buffer(zeros(1))
#     res[1] += x
#     res[1]
# end
# Zygote.gradient(testf, 0.03)
