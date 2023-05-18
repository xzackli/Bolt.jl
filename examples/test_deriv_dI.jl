using DataInterpolations
using ForwardDiff


# elaborate way of writing f(x) = x²
function f(c)
    xs = 0.0:2.0:4.0
    ys = [c * x for x in xs]
    itp = CubicSpline(ys, xs)
    return itp(c)
end


# return DataInterpolations.derivative(itp, 0.5)
    # return itp(0.5)

f(1.0), f(2.0), f(3.0)

##
using DataInterpolations, ForwardDiff

ForwardDiff.derivative(f, 2.0)

##

f_quadratic_spline = c -> square(c, QuadraticSpline)
f_cubic_spline = c -> square(c, CubicSpline)

@test ForwardDiff.derivative(f_quadratic_spline, -3.0) ≈ -6.0
@test ForwardDiff.derivative(f_quadratic_spline, 8.0) ≈ 16.0
@test ForwardDiff.derivative(f_cubic_spline, -3.0) ≈ -6.0
@test ForwardDiff.derivative(f_cubic_spline, 8.0) ≈ 16.0

##

ForwardDiff.derivative(CubicSpline([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]), 0.5)


##
using Test
@test ForwardDiff.derivative(f, -3.0) ≈ -6.0
@test ForwardDiff.derivative(f, 8.0) ≈ 16.0

##
using LinearAlgebra
function g(c)
    B = [c, c + 1]
    A = [
        sin(c)   0.0;
        0.0      cos(c+0.2)
    ]
    x = Tridiagonal(A) \ B
    return x[1]
end

##
g(0.4)


##
ForwardDiff.derivative(g, 0.4)


##

using DataInterpolations, Interpolations, BenchmarkTools
N = 400
u = rand(N)
t = range(0.0, 1.0, length = N)

x = sort(rand(100)) # sorting required for Interpolations

println("cubic spline with DataInterpolations")
interp = DataInterpolations.CubicSpline(u,t)
@profview @btime $interp.($x)
@btime $interp(0.5)
println("cubic spline with Interpolations")
interp2 = Interpolations.CubicSplineInterpolation(t, u)
@btime $interp2.($x)
@btime $interp2(0.5)

println("linear with DataInterpolations")
interplin = DataInterpolations.LinearInterpolation(u,t)
@btime $interplin.($x)
@btime $interplin(0.5)
println("linear with Interpolations")
interplin2 = Interpolations.LinearInterpolation(t, u)
@btime $interplin2.($x)
@btime $interplin2(0.5)
