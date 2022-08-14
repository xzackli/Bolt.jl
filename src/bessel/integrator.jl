

"""
integrates a quadratic polynomial with moments. input f, f′, f″ are those evaluated at a, over
interval (a,b). Uses a preconstructed interpolator for the moments.
"""
function integrate_sph_bessel_filon(f, f′, f″, k, a, b, itp)
    c₂ = f″ / 2
    af″ = a * f″
    c₁ = f′ - af″
    c₀ = f - a * (f′ - af″ / 2)
    k⁻¹ = 1 / k
    k⁻² = k⁻¹ * k⁻¹ 
    k⁻³ = k⁻² * k⁻¹
    
    ΔI = itp(k * b) - itp(k * a)
    s = (c₀ * k⁻¹) * ΔI[1] + (c₁ * k⁻²) * ΔI[2]  + (c₂ * k⁻³) * ΔI[3]
    return s
end

"""
same as above but suitable for use in a piecewise integrator. 
re-use the moment I(ka) which was I(kb) previously
"""
function _loop_integrate_sph_bessel_filon(f, f′, f″, k, a, b, itp, itp_ka)
    c₂ = f″ / 2
    af″ = a * f″
    c₁ = f′ - af″
    c₀ = f - a * (f′ - af″ / 2)
    k⁻¹ = 1 / k
    k⁻² = k⁻¹ * k⁻¹ 
    k⁻³ = k⁻² * k⁻¹
    
    itp_kb = itp(k * b)
    ΔI = itp_kb - itp_ka  # reuse last interpolator
    s = (c₀ * k⁻¹) * ΔI[1] + (c₁ * k⁻²) * ΔI[2]  + (c₂ * k⁻³) * ΔI[3]
    return s, itp_kb
end

