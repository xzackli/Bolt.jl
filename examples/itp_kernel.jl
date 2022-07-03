using CUDA,Interpolations,Adapt

#construct a Bspline
x = 0:0.01:1
y = sin.(x)
spline(f, x_grid) = scale(interpolate(f, BSpline(Cubic(Line(OnGrid())))), x_grid)
s  = spline([yi for yi in y], x)
#^This is the way we do it in the background and ionization history...before saving them as attributes of those objects

#so what actually happens here? He takes the Bspline struct, which normally allows coeffs to be anything you could imagine
#This freedom is a problem for the gpu, so we adapt the bspline into something that doesn't allow the type of the coeffs to be free
#i.e. the Bspline is instantiated i the last line, but 
function Adapt.adapt_structure(to, itp::Interpolations.BSplineInterpolation{T,N,<:Any,IT,Axs}) where {T,N,IT,Axs}
    coefs = Adapt.adapt_structure(to, itp.coefs)
    Tcoefs = typeof(coefs)
    Interpolations.BSplineInterpolation{T,N,Tcoefs,IT,Axs}(coefs, itp.parentaxes, itp.it)
end

#the struct for BSplineInterpolation
# struct BSplineInterpolation{T,N,TCoefs<:AbstractArray,IT<:DimSpec{BSpline},Axs<:Tuple{Vararg{AbstractUnitRange,N}}} <: AbstractInterpolation{T,N,IT}
#     coefs::TCoefs
#     parentaxes::Axs
#     it::IT
# end
#we need to now do this for "scale"
#The struct for scaled Interpolation
# struct ScaledInterpolation{T,N,ITPT,IT,RT} <: AbstractInterpolationWrapper{T,N,ITPT,IT}
#     itp::ITPT
#     ranges::RT
# end

function Adapt.adapt_structure(to, itp::Interpolations.ScaledInterpolation{T,N,ITPT,IT,<:Any}) where {T,N,ITPT,IT}
    s = Adapt.adapt_structure(to,itp.itp)
    Titp = typeof(s)
    ranges = Adapt.adapt_structure(to,itp.ranges)
    RT=typeof(ranges)
    Interpolations.ScaledInterpolation{T,N,Titp,IT,RT}(s,ranges)
end
#Test kernel that just takes the already-called bspline and updates an array
let 
    x′ = π/8
    du = cu([NaN])
    function spl_eval_kernel!(du,s,x)
        du[1]=s(x) # calling the spline inside fails, but passing an interpolated value works (and is analogous to the real example)
        return nothing
    end
    @cuda spl_eval_kernel!(du,cu(s),x′)
    du
    
end
