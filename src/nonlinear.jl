function halofit_params(pk)
    f(x) = √_S(pk, x) - 1
    R = find_zero(f, (1e-2, 1e2), Bisection())
    S_val = _S(pk, R)
    #Although S should be one, I prefer to evaluate it to be consistent with the root
    #finding algo result. S is a callable, while S_val is a fixed value
    dlnSdlnR = R / S_val * _dSdR(pk, R)
    d²lnSdlnR² = dlnSdlnR - dlnSdlnR^2 + (R^2 / S_val) * _d²SdR²(pk, R)
    kσ = 1 / R
    neff = -3 - dlnSdlnR
    C = -d²lnSdlnR²
    return [kσ, neff, C]
end

function _an(neff, C)
    return 10^(1.5222 + 2.8553neff + 2.3706 * neff^2 + 0.9903 * neff^3 + 0.2250 * neff^4 -
    0.6038C) #not added the w₀wₐ term, which is 0.1749ΩDEz(1 + wz)
end

function _bn(neff, C)
    return 10^(-0.5642 + 0.5864neff + 0.5716 * neff^2 - 1.5474C)
    # as above not added the w₀wₐ term, +0.2279*ΩDEz*(1.+w))
end

function _cn(neff, C)
    return 10^(0.3698 + 2.0404neff + 0.8161 * neff^2 + 0.5869C)
end

function _γn(neff, C)
    return 0.1971 - 0.0843neff + 0.8460C
end

function _αn(neff, C)
    return abs(6.0835 + 1.3373neff - 0.1959 * neff^2 - 5.5274C)
end

function _βn(neff, C; fνz=0)
    return 2.0379 - 0.7354neff + 0.3157 * neff^2 + 1.2490 * neff^3 +
    0.3980 * neff^4 - 0.1682C + fνz * (1.081 + 0.395 * neff^2)
end

function _νn(neff)
    return 10^(5.2105 + 3.6902neff)
end

function _n_coefficients(neff, C, fνz)
    an = _an(neff, C)
    bn = _bn(neff, C)
    cn = _cn(neff, C)
    γn = _γn(neff, C)
    αn = _αn(neff, C)
    βn = _βn(neff, C, fνz = fνz)
    μn = 0
    νn = _νn(neff)
    return [an, bn, cn, γn, αn, βn, μn, νn]
end

function _f1a(Ωmz)
    return Ωmz^-0.0732
end

function _f2a(Ωmz)
    return Ωmz^-0.1423
end

function _f3a(Ωmz)
    return Ωmz^0.0725
end

function _f1b(Ωmz)
    return Ωmz^-0.0307
end

function _f2b(Ωmz)
    return Ωmz^-0.0585
end

function _f3b(Ωmz)
    return Ωmz^0.0743
end

function _f_coefficients(Ωmz, frac)
    f1a = _f1a(Ωmz)
    f2a = _f2a(Ωmz)
    f3a = _f3a(Ωmz)
    f1b = _f1b(Ωmz)
    f2b = _f2b(Ωmz)
    f3b = _f3b(Ωmz)
    f1 = frac*f1b+(1-frac)*f1a
    f2 = frac*f2b+(1-frac)*f2a
    f3 = frac*f3b+(1-frac)*f3a
    return [f1, f2, f3]
end

function _Δ²H(an, bn, cn, γn, μn, νn, f1, f2, f3, y)
    return (an * y^(3 * f1)/(1 + bn * y^f2 + (cn * f3 * y)^(3-γn)))/(1 + μn / y + νn / y^2)
    #TODO: μn is fixed to zero. Maybe we can remove it?
end

function _Δ²H(an, bn, cn, γn, μn, νn, f1, f2, f3, y, fνz)
    return _Δ²H(an, bn, cn, γn, μn, νn, f1, f2, f3, y)*(1+fνz*0.977)
end

function _Δ²L(pk, k)
    return k^3 / (2 * π^2) * 10^pk(log10(k))
end

function _Δ²̃L(pk, k, fνz, hz)
    return _Δ²L(pk, k) * (1+fνz*47.48* (k/hz)^2/(1+1.5*(k/hz)^2))
end

function _Δ²Q(Δ²L, αn, βn, y)
    return Δ²L * (1 + Δ²L)^βn / (1 + αn * Δ²L) * exp(-y / 4 - y^2 / 8)
end

function takahashi2012(pk, p::Array{T,1}, Ωmz::Real, k::Real) where {T<:Real}
    #based on Takahashi 2012
    kσ, neff, C = p
    an = _an(neff, C)
    bn = _bn(neff, C)
    cn = _cn(neff, C)
    γn = _γn(neff, C)
    αn = _αn(neff, C)
    βn = _βn(neff, C)
    μn = 0
    νn = _νn(neff)
    f1 = _f1b(Ωmz)
    f2 = _f2b(Ωmz)
    f3 = _f3b(Ωmz)
    y = k / kσ
    Δ²H = _Δ²H(an, bn, cn, γn, μn, νn, f1, f2, f3, y)
    Δ²L = _Δ²L(pk, k)
    Δ²Q = _Δ²Q(Δ²L, αn, βn, y)
    pk_nl = (Δ²Q + Δ²H) * 2 * π^2 / k^3
    return pk_nl
end

function takabird(pk, p::Array{T,1}, Ωmz, Ωrz, Mν, hz, k) where {T<:Real}
    #based on the revised Takahashi2012 present on the Class repository
    Ωνz = Mν/(93.14 * hz^2)
    fνz = Ωνz/Ωmz
    ΩDEz = 1-Ωmz-Ωνz-Ωrz
    #TODO: ask if Bolt Ω_r includes neutrinos or only radiation
    frac = ΩDEz/(1-Ωmz)
    # TODO: check if is required a z dependent fν, Ων etc
    kσ, neff, C = p
    an, bn, cn, γn, αn, βn, μn, νn = _n_coefficients(neff, C, fνz)
    f1, f2, f3 = _f_coefficients(Ωmz, frac)
    y = k / kσ
    Δ²H = _Δ²H(an, bn, cn, γn, μn, νn, f1, f2, f3, y, fνz)
    Δ²̃L = _Δ²̃L(pk, k, fνz, hz)
    Δ²Q = _Δ²Q(Δ²̃L, αn, βn, y)
    pk_nl = (Δ²Q + Δ²H) * 2 * π^2 / k^3
    return pk_nl
end

function takabird(pL, ks, 𝕡)
    InterpPmm = Interpolations.interpolate( log10.(vcat(pL...)),
    BSpline(Cubic(Line(OnGrid()))))
    x = LinRange(log10(first(ks)), log10(last(ks)), length(ks))
    InterpPmm = scale(InterpPmm, x)
    InterpPmm = Interpolations.extrapolate(InterpPmm, Line())
    params_halofit = halofit_params(InterpPmm)
    pNL_takabird = [takabird(InterpPmm, params_halofit, 𝕡.Ω_b+𝕡.Ω_c, 𝕡.Ω_r, 𝕡.Σm_ν, 𝕡.h, k)
                    for k in ks]
    return pNL_takabird
end

function plin_and_nonlin(ks,𝕡,bg,ih,n_q,ℓᵧ,ℓ_ν,ℓ_mν,x)
    pL = [plin(k,𝕡,bg,ih,n_q,ℓᵧ,ℓ_ν,ℓ_mν,x) for k in ks]
    pNL = takabird(pL, ks, 𝕡)
    return vcat(pL...), vcat(pNL...)
end
