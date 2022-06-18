using Interpolations, CUDA, CUDA.Adapt,BenchmarkTools

y = rand(10)
itp = interpolate(y, BSpline(Quadratic(Reflect(OnCell()))))

function Adapt.adapt_structure(to, itp::Interpolations.BSplineInterpolation{T,N,<:Any,IT,Axs}) where {T,N,IT,Axs}
    coefs = Adapt.adapt_structure(to, itp.coefs)
    Tcoefs = typeof(coefs)
    Interpolations.BSplineInterpolation{T,N,Tcoefs,IT,Axs}(coefs, itp.parentaxes, itp.it)
end

function kernel(itp, result)
    result[1] = itp(2.5)
    return nothing
end

result = cu([NaN])

@cuda kernel(cu(itp), result)

result[1] ≈ itp(2.5) # true

y = rand(10)
itp = interpolate(y, BSpline(Quadratic(Reflect(OnCell()))))
cu_itp = cu(itp);
new_x = collect(range(1,10,length=1000000))
cu_new_x = cu(new_x)

@btime itp.(new_x) # 2.5ms
@btime CUDA.@sync cu_itp.(cu_new_x) # 46μs