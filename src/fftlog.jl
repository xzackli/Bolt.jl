
# export plan_fftlog, mul!, fftlogfreq

struct FFTLogPlan{T, OT, AA<:AbstractArray{Complex{T},1}, PT<:Plan, IPT<:Plan}
    L::T
    N::Int
    μ::OT
    q::T
    r₀::T
    k₀r₀::T
    uₘ::AA
    fftplan!::PT
    ifftplan!::IPT
end

function plan_fftlog(r::AA, μ, q, k₀r₀=1.0;
                     kropt=true) where {T, AA<:AbstractArray{T,1}}
    logrmin = log(first(r))
    logrmax = log(last(r))
    r₀ = exp((logrmin + logrmax)/2)
    @assert logrmin < logrmax
    N = length(r)
    L = logrmax - logrmin
    if kropt
        k₀r₀ = k₀r₀_low_ringing(N, μ, q, L, k₀r₀)
    end
    m = fftfreq(N, N)  # get indicies that go from [-N/2] to [N/2]
    uₘ_coeff = similar(r, Complex{T})
    for i in eachindex(m)
        uₘ_coeff[i] = uₘ(m[i], μ, q, L, k₀r₀)
    end
    uₘ_coeff[N÷2+1] = real(uₘ_coeff[N÷2+1])  # eq 19
    fftplan! = plan_fft!(uₘ_coeff)
    ifftplan! = plan_ifft!(uₘ_coeff)
    return FFTLogPlan(L, N, μ, q, r₀, k₀r₀, uₘ_coeff, fftplan!, ifftplan!)
end


U_μ(μ, x) = 2^x * gamma((μ + 1 + x)/2) / gamma((μ + 1 - x)/2)
uₘ(m, μ, q, L, k₀r₀) = (k₀r₀)^(-2π * im * m / L) * U_μ(μ, q + 2π * im * m / L)

function k₀r₀_low_ringing(N, μ, q, L, k₀r₀=1.0)
    # from pyfftlog
    dlnr = L / (N)
    xp = (μ + 1 + q) / 2
    xm = (μ + 1 - q) / 2
    y = π * im / 2 / dlnr
    zp = lgamma(xp + y)
    zm = lgamma(xm + y)
    arg = log(2 / k₀r₀) / dlnr + imag(zp + zm) / π
    return k₀r₀ * exp((arg - round(arg))* dlnr)
end


function mul!(Y, pl::FFTLogPlan, A)
    Y .= A
    pl.fftplan! * Y
    Y .*= pl.uₘ
    pl.ifftplan! * Y
end

function ldiv!(Y, pl::FFTLogPlan, A)
    Y .= A
    pl.fftplan! * Y
    Y ./= pl.uₘ
    pl.ifftplan! * Y
end

function fftlogfreq(pl::FFTLogPlan)
    k₀ = pl.k₀r₀ / pl.r₀
    Nhalf = pl.N ÷ 2
    n = range(-Nhalf,Nhalf,length=pl.N)
    return reverse(k₀ .* exp.(n .* pl.L / pl.N))
end
