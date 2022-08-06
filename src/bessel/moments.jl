# Compute moments x^α J_nu over (0,x)

# Lommel function asymptotic for large t
# we use α - 3/2 for speed when computing the moments xᵏ jᵥ with k∈(0,1,2). pow is slow
function s⁽²⁾(t::T, α_minus_half, ν; kmax=20, tol=1e-16) where T <: Real
    s = one(T)
    sₖ = one(T)
    t⁻¹ = one(T) / t
    t⁻² = t⁻¹^2
    ν² = ν^2
    αm1 = α_minus_half - one(T)/2
    @fastmath for k in 0:kmax
        sₖ *= (ν² - (αm1 - 2k)^2) * t⁻²
        s += sₖ
        if abs(sₖ) < tol * abs(s)
            break
        end
    end
    return s * t^α_minus_half * √(t⁻¹)
end

# Bessel J moment generated from Weniger transformation of hypergeometric ₁F₂
function J_moment_weniger_₁F₂(x::T, ν, α, cache) where T
    return (1/(α+ν+1)) * exp((α+ν+1) * log(x) - ν * log(T(2)) - _lgamma(ν+1)) * 
         weniger1F2((1+α+ν)/2, SA[(3+α+ν)/2, 1+ν], -x^2/4, cache)
end

# Bessel J moment generated from asymptotics of the Lommel function
function J_moment_asymptotic_lommel(x::T, ν, α_minus_half) where T
    α = α_minus_half + one(T)/2
    return J_moment_asymptotic_lommel_prefactor(T, ν, α) + x * (
        (α + ν - 1) * besselj(ν, x) *  s⁽²⁾(x, α_minus_half-1, ν-1) -
        besselj(ν-1, x) *  s⁽²⁾(x, α_minus_half, ν))
end

J_moment_asymptotic_lommel_prefactor(::Type{T}, ν, α) where T = exp(
    log(T(2)) * α + _lgamma((ν + α + 1)/2) - _lgamma((ν - α + 1)/2))

# Bessel J from specialized Lommel function asymptotic for J_nu with ν = 2 + 1/2
@muladd function J_moment_asymptotic_nu_five_halves(x::T, α_minus_half, prefactor) where T
    sinx = sin(x)
    cosx = cos(x)
    x⁻¹ = one(T) / x
    c₁ = sqrt(2/T(π) * x⁻¹)
    besselj_3_2 = sinx * x⁻¹ - cosx
    besselj_5_2 = 3 * sinx * x⁻¹ * x⁻¹ - sinx - 3 * cosx * x⁻¹
    return prefactor + c₁ * x * (
        (α_minus_half + 2) * besselj_5_2 * s⁽²⁾(x, α_minus_half-1, T(3)/2) -
        besselj_3_2 * s⁽²⁾(x, α_minus_half,   T(5)/2))
end


# Spherical Bessel Functions: 
# moments xᵏ jᵥ(x) over (0,x)

@inline sph_j_moment_asymptotic_lommel_prefactor(::Type{T}, ν, k) where T =
    J_moment_asymptotic_lommel_prefactor(T, ν + one(T)/2, k - one(T)/2) * √(T(π)/2)

# spherical Bessel j moment generated from asymptotics of the Lommel function
function sph_j_moment_asymptotic_lommel(x::T, ν, k) where T
    return J_moment_asymptotic_lommel(x, ν + one(T)/2, k - one(T)) * √(T(π)/2)
end

# spherical Bessel J moment generated from Weniger transformation of hypergeometric ₁F₂
function sph_j_moment_weniger_₁F₂(x::T, ν, α, cache) where T
    ν′ = ν + T(3)/2
    return (1/(α+ν+1)) * exp((α+ν+1) * log(x) - (ν′-1) * log(T(2)) - _lgamma(ν′)) * 
         weniger1F2(T(1+α+ν)/2, SVector{2,T}((3+α+ν)/2, (ν′)), -x^2/4, cache) * √(T(π)/2)
end

@muladd function sph_j_moment_asymptotic_nu_five_halves(x::T, k, prefactor) where T
    α_minus_half = k - 1
    sinx = sin(x)
    cosx = cos(x)
    x⁻¹ = one(T) / x
    c₁ = sqrt(x⁻¹)
    besselj_3_2 = sinx * x⁻¹ - cosx
    besselj_5_2 = 3 * sinx * x⁻¹ * x⁻¹ - sinx - 3 * cosx * x⁻¹
    return prefactor + c₁ * x * (
        (α_minus_half + 2) * besselj_5_2 * s⁽²⁾(x, α_minus_half-1, T(3)/2) -
        besselj_3_2 * s⁽²⁾(x, α_minus_half,   T(5)/2))
end
