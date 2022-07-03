using CUDA,Interpolations,Adapt,ForwardDiff

# re-use Bspline code
function Adapt.adapt_structure(to, itp::Interpolations.BSplineInterpolation{T,N,<:Any,IT,Axs}) where {T,N,IT,Axs}
    coefs = Adapt.adapt_structure(to, itp.coefs)
    Tcoefs = typeof(coefs)
    Interpolations.BSplineInterpolation{T,N,Tcoefs,IT,Axs}(coefs, itp.parentaxes, itp.it)
end
function Adapt.adapt_structure(to, itp::Interpolations.ScaledInterpolation{T,N,ITPT,IT,<:Any}) where {T,N,ITPT,IT}
    s = Adapt.adapt_structure(to,itp.itp)
    Titp = typeof(s)
    ranges = Adapt.adapt_structure(to,itp.ranges)
    RT=typeof(ranges)
    Interpolations.ScaledInterpolation{T,N,Titp,IT,RT}(s,ranges)
end

# Our custom interp holder struct
abstract type AbstractHolder{T, U, V, GT} end
struct InterpHolder{T, U, AT<:AbstractArray{U,1}, GT, V, IT<:AbstractInterpolation{V,1}} <: AbstractHolder{T, U, V,GT}
    x::T
    G::GT
    A::AT
    I::IT
end
Adapt.@adapt_structure InterpHolder

# build a the IH struct in a function of the same name (as we would the bg)
function InterpHolder(a::T,xgrid=0:0.01:1) where T
    y = a*sin.(xgrid)
    s  = spline([yi for yi in y], xgrid)
    q = ones(2) #its own (float) type
    ih = InterpHolder(a,xgrid,q,s)
    return ih
end


# extremely simple example with no cuda kernel, just testing forwarddiff array assignment
function t_pbg(x::DT) where DT
    y=[NaN]
    y=convert(Array{DT,1},y) #This line is necessary
    # Otherwise we are assigning a dual type value to a Array{Float32,1}/Array{Float64,1}/CuArray type
    # This is not allowed (type invariance) 
    println("type output: ", typeof(y), " type input:",typeof(x))
    y[1] = x
    return y
end

t_pbg(1.)
t_pbg(ForwardDiff.Dual(1.0))


# Test a simple function that just updates the value of a CuArray
function g_2xdual!(u::AbstractArray{T},x::T) where T
    u[1] = x
    return nothing
end
# A function we want to forward diff through where the kernel is the above simple assignment
function g_pbg(x::DT) where DT
    u = cu([NaN])
    u = cu(convert(AbstractArray{DT,1},u)) # This is necessary for type stability of this function
    println("type check: ", typeof(u) <: AbstractArray{DT}, " and x: ", typeof(x) <: DT) #debug statements for type stability
    println("typeof u ",typeof(u), "eltype of u ",eltype(u), ", typeof x ",typeof(x), ", cux: ",typeof(cu(x)))
    @cuda g_2xdual!(u,x)
    return u
end

# Test that the forward diff works
g_pbg(1.f0) # sanity check that function works
g_pbg(ForwardDiff.Dual(1.f0)) # test on a dual
ForwardDiff.derivative(g_pbg,1.f0) # do the scalar derivative

# Same as above, but now we test on an interpolation call using our InterpHolder struct
spline(f, x_grid) = scale(interpolate(f, BSpline(Cubic(Line(OnGrid())))), x_grid)
function g_interp!(u::AbstractArray{T},ih,x::T) where T
    u[1] = ih.I(x)
    return nothing
end
function ig_pbg(x::DT) where DT
    ih = InterpHolder(x)
    u = cu([NaN])
    u = cu(convert(AbstractArray{DT,1},u))
    @cuda g_interp!(u,cu(ih),x)
    return u
end

# Test the function and derivative
isa(ig_pbg(1.f0), Union{Real,AbstractArray}) # Needs to be true if function is forwarddiffable
ig_pbg(1.f0)
ig_pbg(ForwardDiff.Dual(1.f0))
ForwardDiff.derivative(ig_pbg,1.f0)

# Also independently check the case where we don't require the input point of the spline to be of type T (which will be Dual)
function g2_interp!(u::AbstractArray{T},ih,x) where T
    u[1] = ih.I(x)
    return nothing
end
function ig2_pbg(θ::DT) where DT
    ih = InterpHolder(θ)
    u = cu([NaN])
    u = cu(convert(AbstractArray{DT,1},u))
    x = 1.f0
    @cuda g2_interp!(u,cu(ih),x)
    return u
end
ig2_pbg(1.f0)
ig2_pbg(ForwardDiff.Dual(1.f0))
ForwardDiff.derivative(ig2_pbg,1.f0)
# And it does still work - u is the only thing that impacts type stability (which makes sense)


println("Done.")
# Old stuff ----

# What we actally want here is the following.

# f(𝕡) = plin_𝕡(𝕡)
# ∂𝕡f(𝕡2) = ForwardDiff.derivative(f,𝕡2)
# But we have a function that looks like 
# function plin_𝕡(𝕡)
    # q = g(𝕡)
    # ...
    # some stuff on the cpu...
    # ...
    # @cuda gpu_y(...g(𝕡), z)
    # ...
    # a = h(z)
    # return a
# end

# So we *don't* want to take forward diffs INSIDE a gpu kernel (at least, not here)
# Instead we want to take forward diffs THROUGH a gpu kernel (which is probably easier?)
# So let's make a toy example of a forward diff of a function that contains a gpu kernel
# and check that the ForwardDiff just works through it. Then move up to interpolations and structs.
# Basically we just do what we've already been doing but wrap everything with a function that gets 
# ForwardDiff'd.
using Cthulhu

function gpu_sin!(y1,x1,a1)
    y1[1] = a1*x1
    return nothing
end
gpu_sin!([NaN],2,3)

let 
function gpu_sin!(y1,x1,a1)
    y1[1] = a1*sin(x1)
    return nothing
end

function f(a2::T) where T
    x2 =  1.
    y2 = cu([NaN])
    #@device_code_warntype interactive=true 
    @cuda gpu_sin(y2,T(x2),T(a2))
    return y2
end
println(f(ForwardDiff.Dual(1.)))
# ForwardDiff.derivative(f, 1.)
end