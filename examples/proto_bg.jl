using CUDA,Interpolations,Adapt

#re-use Bspline code
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

#this almost works...fails on adapt I THINK because AT cannot be changed after the fact from the inferred (Float64) array to a CuArray by adapt?
#why does this hppae
# abstract type AbstractArrayHolder{T} end
# struct ArrayHolder{T} <: AbstractArrayHolder{T} 
#     A::AbstractArray{T,1}
# end
abstract type AbstractArrayHolder{T} end
struct ArrayHolder{T,AT<:AbstractArray{T,1}} <: AbstractArrayHolder{T} 
    A::AT
end
function Adapt.adapt_structure(to,ah::ArrayHolder{T,AT}) where {T,AT}
    A=Adapt.adapt_storage(to,ah.A)
    t=eltype(A)
    at=typeof(A) 
    ArrayHolder{t,at}(A)
end

# Need the extra type V here, there is no reason the interpolation and array both share types
# And the example below breaks if I don't do this...
abstract type AbstractInterpArrayHolder{T,V} end
struct InterpArrayHolder{T,AT<:AbstractArray{T,1},V,IT<:AbstractInterpolation{V,1}} <: AbstractInterpArrayHolder{T,V} 
    A::AT
    I::IT
end
function Adapt.adapt_structure(to,ia::InterpArrayHolder{T,AT,V,IT}) where {T,AT,V,IT}
    A=Adapt.adapt_storage(to,ia.A)
    at=typeof(A) 
    I=Adapt.adapt_structure(to,ia.I)
    it=typeof(I)
    t=eltype(A)
    v=eltype(I)
    # InterpArrayHolder{t,at,v,it}(A,I)
    InterpArrayHolder(A,I)
end


abstract type AbstractHolder{T, U, V, GT} end
struct InterpHolder{T, U, AT<:AbstractArray{U,1}, GT, V, IT<:AbstractInterpolation{V,1}} <: AbstractHolder{T, U, V,GT}
    x::T
    G::GT
    A::AT
    I::IT
end


function Adapt.adapt_structure(to, ih::InterpHolder{T, U, AT, GT, V,IT}) where {T, U, AT, GT, V,IT}
    #adapt the fields
    I = Adapt.adapt_structure(to,ih.I)
    it = typeof(I)
    G = Adapt.adapt_structure(to,ih.G)
    gt = typeof(G)
    x = Adapt.adapt_structure(to,ih.x)
    A = Adapt.adapt_storage(to,ih.A)
    at=typeof(A)
    t=typeof(x)
    u=eltype(A)
    v=eltype(I)
    # InterpHolder{t,u,at,gt,v,it}(x,G,A,I)
    InterpHolder(x,G,A,I)
end

#construct a Bspline and store it in our interp_holder struct
x = 0:0.01:1
y = sin.(x)
spline(f, x_grid) = scale(interpolate(f, BSpline(Cubic(Line(OnGrid())))), x_grid)
s  = spline([yi for yi in y], x)
q = ones(2)
ih = InterpHolder(0.,x,q,s)
ia = InterpArrayHolder(q,s)
ih.I(.5)
typeof(q) <: AbstractArray#{Float64,1}
# ah = ArrayHolder{eltype(q),typeof(q)}(q) #this works
CuArray<:AbstractArray
ah = ArrayHolder(q)
ah.A
eltype(s)

#test macro
Adapt.@adapt_structure InterpHolder

#Test kernel that queries the scaled b-spline, which is an attribute of an "interp_holder" struct
let    
    x′ = π/8
    du = cu([NaN])
    function ih_spl_eval_kernel!(du,ih,x)
    # function ia_spl_eval_kernel!(du,ia,x)
    # function a_kernel!(du,ah,x)
        du[1]=ih.I(x)
        # du[1]=ia.I(x) 
        # du[1]=ia.A[1] 
        # du[1]=ah.A[1]
        return nothing
    end
    @cuda ih_spl_eval_kernel!(du,cu(ih),x′)
    # @cuda ia_spl_eval_kernel!(du,cu(ia),x′)
    # @cuda a_kernel!(du,cu(ah),x′)
    du
    
end
