# Generate tables to interpolate for the moments of spherical bessel functions
# This is mostly to avoid ₁F₂ for kη < 50, but also still provides a nice 7x speedup over 
# evaluating the asymptotic derived from the Lommel function.
# TODO: use maclaurin series for small arguments?
# TODO: store interpolation tables as artifacts?

# parametrize the internal type T, the output type TX, and the Bessel multipole NU
struct MomentTableConfig{T, TX, NU, ORDER}
    ν::Int
end


sph_j_moment_asymp(config::MomentTableConfig{T, TX, 2}, x, m, pref) where {T,TX} =
    sph_j_moment_asymp_nu_2(x, m, pref)
sph_j_moment_asymp(config::MomentTableConfig{T, TX, 3}, x, m, pref) where {T,TX} =
    sph_j_moment_asymp_nu_3(x, m, pref)
sph_j_moment_asymp(config::MomentTableConfig, x, m, pref) = sph_j_moment_asymp(x, config.ν, m, pref)

# it's important for this general function to be fast, for source function integration.
# For general Cₗ, we can't afford to store a full interpolator.
# TODO: we can make this 2x faster by removing redundant calculations between the different m
@inline function sph_j_moment_asymp_third_order(config::MomentTableConfig{T, TX, V, 3}, x,
                                          pref0, pref1, pref2) where {T,TX,V}
    return SVector{3,TX}(
        sph_j_moment_asymp(config, (x), 0, pref0),
        sph_j_moment_asymp(config, (x), 1, pref1),
        sph_j_moment_asymp(config, (x), 2, pref2))
end

@inline function sph_j_moment_weniger_third_order(config::MomentTableConfig{T, TX, V, 3}, x, cache) where {T,TX,V}
    ν = config.ν
    return SVector{3,TX}(
        sph_j_moment_weniger_₁F₂(T(x), ν, 0, cache),
        sph_j_moment_weniger_₁F₂(T(x), ν, 1, cache),
        sph_j_moment_weniger_₁F₂(T(x), ν, 2, cache))
end

# TODO: fourth order moments, i.e. x^3 terms

struct MomentTable{T, TX, ITP, C, W}
    itp::ITP
    config::C
    pref0::T
    pref1::T
    pref2::T
    kη_min::T
    kη_max::T
    cache::W
    c₁::TX
    c₂::TX
end

function (m::MomentTable{T,TX})(x::TX) where {T,TX}
    if m.kη_min ≤ x ≤ m.kη_max
        return m.itp(x)::SVector{3,TX}
    elseif x > m.kη_max
        return sph_j_moment_asymp_third_order(
            m.config, x, m.pref0, m.pref1, m.pref2)::SVector{3,TX}
    else 
        # return sph_j_moment_weniger_third_order(m.config, x, m.cache)::SVector{3,TX}
        error("Do not extrapolate the moment table downwards. kη_min = $(m.kη_min) > $(x).")
    end
end

Base.show(io::IO, itp::MomentTable{T, TX, ITP, MomentTableConfig{T, TX, NU}}) where {T, TX, ITP, NU} = 
    print(io, "MomentTable($T internal, $TX output, ν=$NU) with kη over $(first(itp.itp.ranges))")

# T is internal precision, TX is output precision
function sph_bessel_interpolator(::Type{T}, ::Type{TX}, ν::Int, max_order, kη_min, kη_max, N::Int; 
                                        weniger_cut=50) where {T, TX}

    pref0 = sph_j_moment_asymp_prefactor(T, ν, 0)
    pref1 = sph_j_moment_asymp_prefactor(T, ν, 1)
    pref2 = sph_j_moment_asymp_prefactor(T, ν, 2)
    c₁ = (1/(α+ν+1) )* √(T(π)/2)
    c₂ = (ν + T(1)/2) * log(T(2)) + _lgamma(ν + T(3)/2)
    cache = WenigerCache1F2(T)

    config = MomentTableConfig{T, TX, ν, max_order}(ν)

    xs = LinRange(kη_min, kη_max, N)
    xs_Tx = LinRange(TX(kη_min), TX(kη_max), N)
    moments = Vector{SVector{3,TX}}(undef, N)
    for i in eachindex(xs)
        x = xs[i]
        if x < weniger_cut
            moments[i] = sph_j_moment_weniger_third_order(config, x, cache)
        else
            moments[i] = sph_j_moment_asymp_third_order(config, x, pref0, pref1, pref2)
        end
    end
    itp = scale(interpolate(moments, BSpline(Cubic(Line(OnGrid())))), xs_Tx)
    return MomentTable{T, TX, typeof(itp), typeof(config), WenigerCache1F2}(
        itp, config, pref0, pref1, pref2, kη_min, kη_max, cache, c₁, c₂)
end
sph_bessel_interpolator(ν::Int, args...; kwargs...) =
    sph_bessel_interpolator(Double64, Float64, ν, args...; kwargs...)
