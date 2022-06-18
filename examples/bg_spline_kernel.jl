using CUDA,ForwardDiff, Adapt, Interpolations
using Pkg
Pkg.activate("/pscratch/sd/j/jsull/julia/Bolt.jl")
using Bolt
using BenchmarkTools

# Start with a sanity check from the runtests.jl file
function fbg(Ω_b::DT) where DT
    𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
    return bg.η(-5.)
end
fbg(0.046)
ForwardDiff.derivative(fbg,0.046)

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


# ----

# The overview is that calling the bg from the GPU currently does not work because we have not adapted
# the bg struct 

# First we need to redefine the background struct to not use "Array" since this doesn't work with CuArrays
# This is done in the background.jl file

# Then we need to make sure that this struct is adapted to be amenable to passing to GPU
Adapt.@adapt_structure Background

# We don't actually need to redefine the Background function with the constructor in it, since it is type stable.


function g_interp_η!(u::AbstractArray{T},bg,x) where T
    u[1] = bg.η(x)
    return nothing
end

function g_fbg(Ω_b::DT) where DT
    𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
    u = cu([NaN])
    u = cu(convert(AbstractArray{DT,1},u))
    x = -5.
    @cuda g_interp_η!(u,cu(bg),x)
    return u
end

# Sanity test that bg still does the normal thing with normal numbers and Duals/ForwardDiff
# Check that the bg can be passsed to a GPU kernel and have its interpolations called from GPU
g_fbg(0.046)
g_fbg(ForwardDiff.Dual(0.046)) # test on a dual

# Compute ForwardDiff derivative on the GPU version
ForwardDiff.derivative(g_fbg,0.046)
