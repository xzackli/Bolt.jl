
using CUDA, CUDA.Adapt
import CUDA.Adapt: adapt_structure

function adapt_structure(to, x::ComponentArray)
    ComponentArray{getaxes(x)}(adapt_structure(to, getdata(x)))
end
function adapt_structure(to, itp::Interpolations.BSplineInterpolation{T,N,<:Any,IT,Axs}) where {T,N,IT,Axs}
    coefs = adapt_structure(to, itp.coefs)
    Interpolations.BSplineInterpolation{T,N,typeof(coefs),IT,Axs}(coefs, itp.parentaxes, itp.it)
end
function adapt_structure(to, itp::Interpolations.ScaledInterpolation{T,N,<:Any,IT,<:Any}) where {T,N,IT}
    s = adapt_structure(to, itp.itp)
    ranges = adapt_structure(to, itp.ranges)
    Interpolations.ScaledInterpolation{T,N,typeof(s),IT,typeof(ranges)}(s, ranges)
end
Adapt.@adapt_structure Background
Adapt.@adapt_structure IonizationHistory
Adapt.@adapt_structure Hierarchy

function CUDA.GPUArrays._copyto!(dest::ComponentArray, bc::Base.Broadcast.Broadcasted)
    CUDA.GPUArrays._copyto!(getdata(dest), bc)
    dest
end
