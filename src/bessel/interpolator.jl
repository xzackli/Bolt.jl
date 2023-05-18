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

# it's important for this general function to be fast, for source function integration.
# For general Cₗ, we can't afford to store a full interpolator.
# TODO: we can make this 2x faster by removing redundant calculations between the different m

const _sjm = sph_j_moment_maclaurin_₁F₂
const _sjw = sph_j_moment_weniger_₁F₂
const _sja = sph_j_moment_asymp

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

struct MomentTable{T, TX, O, C, ITP, W}
    itp::ITP
    config::C
    prefactors::NTuple{O,T}
    kη_min::T
    kη_max::T
    cache::W
    c₂::TX  # (ν + T(1)/2) * log(T(2)) + _lgamma(ν + T(3)/2)
end

# trait to extract ν from type of MomentTable, it's a little buried...
getnu(::MomentTable{T, TX, O, C}) where {T,TX,O,V,C<:MomentTableConfig{T, TX, V}} = V
getorder(::MomentTable{T, TX, O}) where {T,TX,O} = O


# spherical Bessel J moment generated from Maclaurin series of hypergeometric ₁F₂
function sph_j_moment_maclaurin_₁F₂(x::T, m, mt::MomentTable) where T
    ν = getnu(mt)
    return (√(T(π)/2)/(m+ν+1)) * exp((m+ν+1) * log(x) - mt.c₂) * 
        maclaurin_₁F₂(T(1+m+ν)/2, SVector{2,T}((3+m+ν)/2, (ν + T(3)/2)), -x^2/4)
end

# this converts to internal type T (typically Double64) and then converts back to output type Tx)
@inline function sph_j_moment_maclaurin_all_orders(mt::MomentTable{T, TX, O, MomentTableConfig{T, TX, V, 3}}, 
    x) where {T,TX,O,V}
    return SVector{3,TX}(_sjm(TX(x), 0, mt), _sjm(TX(x), 1, mt), _sjm(TX(x), 2, mt))
end

function (mt::MomentTable{T,TX,O})(x::TX) where {T,TX,O}
    if mt.kη_min ≤ x ≤ mt.kη_max
        return mt.itp(x)::SVector{O,TX}
    elseif x > mt.kη_max
        return sph_j_moment_asymp_all_orders(mt.config, x, mt.prefactors...)::SVector{O,TX}
    else 
        return sph_j_moment_maclaurin_all_orders(mt, x)::SVector{O,TX}
    end
end

Base.show(io::IO, itp::MomentTable{T, TX, O}) where {T, TX, O} = 
    print(io, "MomentTable($T internal, $TX output, ν=$(getnu(itp))) of order $O with kη over $(first(itp.itp.ranges))")

# T is internal precision, TX is output precision
function sph_bessel_interpolator(::Type{T}, ::Type{TX}, ν::Int, order, kη_min, kη_max, N::Int; 
                                        weniger_cut=50) where {T, TX}

    prefactors = Tuple(sph_j_moment_asymp_prefactor(T, ν, m) for m in 0:(order-1))
    c₂ = (ν + T(1)/2) * log(T(2)) + _lgamma(ν + T(3)/2)
    cache = WenigerCache1F2(T)
    config = MomentTableConfig{T, TX, ν, order}(ν)

    itp = _fill_sph_bessel_interp(config, kη_min, kη_max, N, cache, weniger_cut, prefactors)
    return MomentTable{T, TX, order, typeof(config), typeof(itp), WenigerCache1F2}(
        itp, config, prefactors, kη_min, kη_max, cache, c₂)
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
    return Interpolations.scale(Interpolations.interpolate(
        moments, Interpolations.BSpline(Interpolations.Cubic(
            Interpolations.Line(Interpolations.OnGrid())))), xs_Tx)
end




## fourth order versions: up to cubic terms

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



@inline function sph_j_moment_maclaurin_all_orders(mt::MomentTable{T, TX, O, MomentTableConfig{T, TX, V, 4}}, 
    x) where {T,TX,O,V}
    return SVector{4,TX}(_sjm(T(x), 0, mt), _sjm(T(x), 1, mt), _sjm(T(x), 2, mt), _sjm(T(x), 3, mt))
end


# we generally always use an interpolator for weniger transformation, so these are un-used

# spherical Bessel J moment generated from Weniger transformation of hypergeometric ₁F₂
# this version stores the constant c₂ = (ν + 1/2) log(2) + Γ(ν + 3/2)
# function sph_j_moment_weniger_₁F₂(x::T, m, mt::MomentTable) where T
#     ν = getnu(mt)
#     return (√(T(π)/2)/(m+ν+1)) * exp((m+ν+1) * log(x) - mt.c₂) * 
#          weniger1F2(T(1+m+ν)/2, SVector{2,T}((3+m+ν)/2, (ν + T(3)/2)), -x^2/4, mt.cache)
# end

# @inline function sph_j_moment_weniger_all_orders(mt::MomentTable{T, TX, O, MomentTableConfig{T, TX, V, 4}}, 
#     x) where {T,TX,O,V}
# return SVector{4,TX}(_sjw(T(x), 0, mt), _sjw(T(x), 1, mt), _sjw(T(x), 2, mt), _sjw(T(x), 3, mt))
# end

# this converts to internal type T (typically Double64) and then converts back to output type Tx)
# @inline function sph_j_moment_weniger_all_orders(mt::MomentTable{T, TX, O, MomentTableConfig{T, TX, V, 3}}, 
#                                                  x) where {T,TX,O,V}
#     return SVector{3,TX}(_sjw(T(x), 0, mt), _sjw(T(x), 1, mt), _sjw(T(x), 2, mt))
# end

# sph_j_moment_asymp(config::MomentTableConfig, x, m, pref) = sph_j_moment_asymp(x, config.ν, m, pref)
