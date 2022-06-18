using CUDA,ForwardDiff, Adapt, Interpolations
using Pkg
Pkg.activate("/pscratch/sd/j/jsull/julia/Bolt.jl")
using Bolt

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


Adapt.@adapt_structure IonizationHistory
function g_interp_csb²!(u::AbstractArray{T},ih,x) where T
    u[1] = ih.csb²(x)
    return nothing
end

function g_fih(Ω_b::DT) where DT
    𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    u = cu([NaN])
    u = cu(convert(AbstractArray{DT,1},u))
    x = -5.
    @cuda g_interp_csb²!(u,cu(ih),x)
    return u
end


g_fih(0.046)
#g_fih(ForwardDiff.Dual(0.046)) # test on a dual - this fails for some reason...
ForwardDiff.derivative(g_fih,0.046)
