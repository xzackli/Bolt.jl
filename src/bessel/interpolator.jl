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

const _sja = sph_j_moment_asymp
const _sjw = sph_j_moment_weniger_₁F₂

# 3rd order: up to quadratic terms
@inline function sph_j_moment_asymp_all_orders(config::MomentTableConfig{T, TX, V, 3}, x,
                                                p0, p1, p2) where {T,TX,V}
    return SVector{3,TX}(_sja(config, x, 0, p0), _sja(config, x, 1, p1), _sja(config, x, 2, p2))
end

# this converts to internal type T (typically Doubble64) and then converts back to output type Tx)
@inline function sph_j_moment_weniger_all_orders(config::MomentTableConfig{T, TX, V, 3}, 
                                                  x, cache) where {T,TX,V}
    ν = config.ν
    return SVector{3,TX}(_sjw(T(x), ν, 0, cache), _sjw(T(x), ν, 1, cache), _sjw(T(x), ν, 2, cache))
end

# 4th order: up to cubic terms
@inline function sph_j_moment_asymp_all_orders(config::MomentTableConfig{T, TX, V, 4}, x,
                                              p0, p1, p2, p3) where {T,TX,V}
    return SVector{4,TX}(_sja(config, x, 0, p0), _sja(config, x, 1, p1), 
        _sja(config, x, 2, p2), _sja(config, x, 3, p3))
end

@inline function sph_j_moment_weniger_all_orders(config::MomentTableConfig{T, TX, V, 4}, 
                                                 x, cache) where {T,TX,V}
    ν = config.ν
    return SVector{4,TX}(_sjw(T(x), ν, 0, cache), _sjw(T(x), ν, 1, cache), 
        _sjw(T(x), ν, 2, cache), _sjw(T(x), ν, 3, cache))
end



struct MomentTable{T, TX, O, C, ITP, W}
    itp::ITP
    config::C
    prefactors::NTuple{O,T}
    kη_min::T
    kη_max::T
    cache::W
    c₁::TX
    c₂::TX
end

function (m::MomentTable{T,TX,O})(x::TX) where {T,TX,O}
    if m.kη_min ≤ x ≤ m.kη_max
        return m.itp(x)::SVector{O,TX}
    elseif x > m.kη_max
        return sph_j_moment_asymp_all_orders(m.config, x, m.prefactors...)::SVector{O,TX}
    else 
        # return sph_j_moment_weniger_all_orders(m.config, x, m.cache)::SVector{3,TX}
        error("Do not extrapolate the moment table downwards. kη_min = $(m.kη_min) > $(x).")
    end
end

Base.show(io::IO, itp::MomentTable{T, TX, O, MomentTableConfig{T, TX, NU, O}}) where {T, TX, O, NU} = 
    print(io, "MomentTable($T internal, $TX output, ν=$NU) of order $O with kη over $(first(itp.itp.ranges))")

# T is internal precision, TX is output precision
function sph_bessel_interpolator(::Type{T}, ::Type{TX}, ν::Int, order, kη_min, kη_max, N::Int; 
                                        weniger_cut=50) where {T, TX}

    prefactors = Tuple(sph_j_moment_asymp_prefactor(T, ν, m) for m in 0:(order-1))
    c₁ = (1/(α+ν+1) )* √(T(π)/2)
    c₂ = (ν + T(1)/2) * log(T(2)) + _lgamma(ν + T(3)/2)
    cache = WenigerCache1F2(T)
    config = MomentTableConfig{T, TX, ν, order}(ν)

    itp = _fill_sph_bessel_interp(config, kη_min, kη_max, N, cache, weniger_cut, prefactors)
    return MomentTable{T, TX, order, typeof(config), typeof(itp), WenigerCache1F2}(
        itp, config, prefactors, kη_min, kη_max, cache, c₁, c₂)
end
sph_bessel_interpolator(ν::Int, args...; kwargs...) =
    sph_bessel_interpolator(Double64, Float64, ν, args...; kwargs...)

function _fill_sph_bessel_interp(config::MomentTableConfig{T, TX, NU, O}, 
                                 kη_min, kη_max, N, cache, weniger_cut, prefactors) where {T, TX, NU, O}
    xs = LinRange(kη_min, kη_max, N)
    xs_Tx = LinRange(TX(kη_min), TX(kη_max), N)
    moments = Vector{SVector{O,TX}}(undef, N)
    for i in eachindex(xs)
        x = xs[i]
        if x < weniger_cut
            moments[i] = sph_j_moment_weniger_all_orders(config, x, cache)
        else
            moments[i] = sph_j_moment_asymp_all_orders(config, x, prefactors...)
        end
    end
    return scale(interpolate(moments, BSpline(Cubic(Line(OnGrid())))), xs_Tx)
end

