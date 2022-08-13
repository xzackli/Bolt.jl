

function integrate_sph_bessel_filon(f, f′, f″, k, a, b, itp)
    c₂ = f″ / 2
    af″ = a * f″
    c₁ = f′ - af″
    c₀ = f - a * (f′ - (a/2) * af″)
    k⁻¹ = 1 / k
    k⁻² = k⁻¹ * k⁻¹ 
    k⁻³ = k⁻² * k⁻¹
    
    ΔI = itp(k * b) - itp(k * a)  # reuse last interpolator
    s = (c₀ * k⁻¹) * ΔI[1] + (c₁ * k⁻²) * ΔI[2]  + (c₂ * k⁻³) * ΔI[3]
    return s
end


function _loop_integrate_sph_bessel_filon(f, f′, f″, k, a, b, itp, itp_ka)
    c₂ = f″ / 2
    af″ = a * f″
    c₁ = f′ - af″
    c₀ = f - a * (f′ - (a/2) * af″)
    k⁻¹ = 1 / k
    k⁻² = k⁻¹ * k⁻¹ 
    k⁻³ = k⁻² * k⁻¹
    
    itp_kb = itp(k * b)
    ΔI = itp_kb - itp_ka  # reuse last interpolator
    s = (c₀ * k⁻¹) * ΔI[1] + (c₁ * k⁻²) * ΔI[2]  + (c₂ * k⁻³) * ΔI[3]
    return s, itp_kb
end

