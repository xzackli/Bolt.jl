
# utilities for x ↔ scale factor ↔ redshift
a2z(a::T) where T = one(T)/a - one(T)
z2a(z::T) where T = one(T)/(z + one(T))
a2x(a) = log(a)
x2a(x) = exp(x)
z2x(z) = a2x(z2a(z))
x2z(x) = a2z(x2a(x))

# utility function for constructing an interpolator
spline(f, x_grid) = scale(interpolate(f, BSpline(Cubic(Line(OnGrid())))), x_grid)
spline_∂ₓ(f, x_grid) = spline([Interpolations.gradient(f, x)[1] for x in x_grid], x_grid)
spline_∂ₓ²(f, x_grid) = spline([Interpolations.hessian(f, x)[1] for x in x_grid], x_grid)

δ_kron(i, j) = (i == j) ? 1 : 0

"""log of the gamma function"""
_lgamma(x) = @inbounds logabsgamma(x)[1]
_lgamma(x::Double64) = lgamma(x)



#utilities for mapping between comoving momenta and unit interval
to_ui(lq,lqmi,lqma) = -1 + (1- (-1)) / (lqma-lqmi) * (lq-lqmi)
from_ui(x,lqmi,lqma) = lqmi + (lqma- lqmi) / (1- (-1)) * (x- (-1))
dxdq(q,logqmin,logqmax) = (1+to_ui(1+logqmin,logqmin,logqmax))/(q*log(10))
xq2q(x,logqmin,logqmax) = 10.0 ^ from_ui(x,logqmin,logqmax)

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))



# utilities for applying an FFTLog transformation ===

struct FFTLogPlan{T, OT, AA<:AbstractArray{Complex{T},1},
                  AAR<:AbstractArray{T,1}, PT<:Plan, IPT<:Plan}
    L::T
    N::Int
    μ::OT
    q::T
    r₀::T
    k₀r₀::T
    uₘ::AA
    r::AAR
    k::AAR
    fftplan!::PT
    ifftplan!::IPT
end

function plan_fftlog(r::AA, μ, q, uₘ, k₀r₀=1.0;
                     kropt=true) where {T, AA<:AbstractArray{T,1}}
    logrmin = log(first(r))
    logrmax = log(last(r))
    r₀ = exp((logrmin + logrmax)/2)
    @assert logrmin < logrmax
    N = length(r)
    L = logrmax - logrmin
    dlnr = L / (N - 1)
    if kropt
        k₀r₀ = k₀r₀_low_ringing(N, μ, q, L, k₀r₀)
    end
    k₀ = k₀r₀ / r₀
    Nhalf = N ÷ 2
    n = range(-Nhalf,Nhalf,length=N)
    k = reverse(k₀ .* exp.(n .* L / N))

    m = fftfreq(N, N)  # get indicies that go from [-N/2] to [N/2]
    uₘ_coeff = similar(r, Complex{T})
    for i in eachindex(m)
        uₘ_coeff[i] = uₘ(m[i], μ, q, dlnr, k₀r₀, N)
    end
    uₘ_coeff[N÷2+1] = real(uₘ_coeff[N÷2+1])  # eq 19
    fftplan! = plan_fft!(uₘ_coeff)
    ifftplan! = plan_ifft!(uₘ_coeff)
    return FFTLogPlan(L, N, μ, q, r₀, k₀r₀, uₘ_coeff, r, k, fftplan!, ifftplan!)
end


U_μ(μ, x) = exp(x * log(2.) - loggamma(0.5 * (μ + 1 - x)) + loggamma(0.5 * (μ + 1 + x)))
uₘ(m, μ, q, dlnr, k₀r₀, N) = (k₀r₀)^(-2π * im * m / (dlnr * N)) * U_μ(μ, q + 2π * im * m / (dlnr * N) )

#3D Bessel j Mellin kernel
M_μ(μ, x) = exp( (x-0.5)*log(2.) - lgamma(0.5 * (μ + 2 - x)) + lgamma(0.5 * (μ + 1 + x)))
mₘ(m, μ, q, dlnr, k₀r₀, N) = (k₀r₀)^(-2π * im * m / (dlnr * N)) * M_μ(μ, q + 2π * im * m / (dlnr * N) )


function k₀r₀_low_ringing(N, μ, q, L, k₀r₀=1.0)
    # from pyfftlog
    dlnr = L / (N-1)
    xp = (μ + 1 + q) / 2
    xm = (μ + 1 - q) / 2
    y = π * im / 2 / dlnr
    zp = loggamma(xp + y)
    zm = loggamma(xm + y)
    arg = log(2 / k₀r₀) / dlnr + imag(zp + zm) / π
    return k₀r₀ * exp((arg - round(arg))* dlnr)
end

function mul!(Y, pl::FFTLogPlan, A)
    Y .= A
    Y .*= (pl.r).^(-pl.q)
    pl.fftplan! * Y
    Y .*= pl.uₘ
    pl.ifftplan! * Y
    Y .*= (pl.r).^(pl.q)
end

function ldiv!(Y, pl::FFTLogPlan, A)
    Y .= A
    Y .*= (pl.r).^(-pl.q)
    pl.fftplan! * Y
    Y ./= pl.uₘ
    pl.ifftplan! * Y
    Y .*= (pl.r).^(pl.q)
end
