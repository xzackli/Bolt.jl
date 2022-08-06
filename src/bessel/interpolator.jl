

## moment Interpolations


@inline function sph_j_moment_weniger_0_1_2(::Type{T}, ::Type{TX}, x, ν, cache) where {T,TX}
    return SVector{3,TX}(
        sph_j_moment_weniger_₁F₂(T(x), ν, 0, cache),
        sph_j_moment_weniger_₁F₂(T(x), ν, 1, cache),
        sph_j_moment_weniger_₁F₂(T(x), ν, 2, cache))
end

@inline function sph_j_moment_asymptotic_nu_2_0_1_2(::Type{T}, ::Type{TX}, x, 
                                                    pref0, pref1, pref2) where {T,TX}
    return SVector{3,TX}(
        sph_j_moment_asymptotic_nu_five_halves(T(x), 0, pref0),
        sph_j_moment_asymptotic_nu_five_halves(T(x), 1, pref1),
        sph_j_moment_asymptotic_nu_five_halves(T(x), 2, pref2))
end


struct MomentTable{T, TX, NU, ITP, W}
    itp::ITP
    pref0::T
    pref1::T
    pref2::T
    kη_min::T
    kη_max::T
    cache::W
    c₁::TX
    c₂::TX
end

function (m::MomentTable{T, TX, 2})(x::TX) where {T, TX}
    if m.kη_min < x < m.kη_max
        return m.itp(x)::SVector{3,TX}
    elseif x > m.kη_max
        return sph_j_moment_asymptotic_nu_2_0_1_2(TX, TX, x, m.pref0, 
                                                  m.pref1, m.pref2)::SVector{3,TX}
    else
        return sph_j_moment_weniger_0_1_2(T, TX, x, 2, m.cache)::SVector{3,TX}
    end
end

Base.show(io::IO, itp::MomentTable{T, TX, NU, ITP, W}) where {T, TX, NU, ITP, W} = 
    print(io, "MomentTable{$T, $TX, $NU} with kη over $(first(itp.itp.ranges))")

# function sph_bessel_dipole_interpolator(::Type{T}, ::Type{TX}, kη_min, kη_max, N::Int; 
#                                         weniger_cut=50) where {T, TX}

#     ν = 2
#     pref0 = sph_j_moment_asymptotic_prefactor(T, ν, 0)
#     pref1 = sph_j_moment_asymptotic_prefactor(T, ν, 1)
#     pref2 = sph_j_moment_asymptotic_prefactor(T, ν, 2)
#     c₁ = (1/(α+ν+1) )* √(T(π)/2)
#     c₂ =  (ν + T(1)/2) * log(T(2)) + _lgamma(ν + T(3)/2)
#     cache = WenigerCache1F2(T)

#     xs = LinRange(kη_min, kη_max, N)
#     xs_Tx = LinRange(TX(kη_min), TX(kη_max), N)
#     moments = Vector{SVector{3,TX}}(undef, N)
#     for i in eachindex(xs)
#         x = xs[i]
#         if x < weniger_cut
#             moments[i] =  sph_j_moment_weniger_0_1_2(T, TX, x, ν, cache)
#         else
#             moments[i] = sph_j_moment_asymptotic_nu_2_0_1_2(T, TX, x, pref0, pref1, pref2)
#         end
#     end
#     itp = scale(interpolate(moments, BSpline(Cubic(Line(OnGrid())))), xs_Tx)
#     return MomentTable{T, TX, 2, typeof(itp), WenigerCache1F2}(
#         itp, pref0, pref1, pref2, kη_min, kη_max, cache, c₁, c₂)

# end
# sph_bessel_dipole_interpolator(kη_min::T, args...; kwargs...) where {T <: Real} =
#     sph_bessel_dipole_interpolator(Double64, Float64, kη_min, args...; kwargs...)

# sph_bessel_dipole_interpolator(::Type{T}, args...; kwargs...) where {T<:Float64} = 
#     error("This won't work with precision less than Double64!")

