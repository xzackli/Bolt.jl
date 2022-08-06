# Adapted from HypergeometricFunctions, this computes ₁F₂ with a sequence transformation 
# from Weniger. I've modified it to avoid allocations so when we generate the moment 
# table, so we can multithread when generating the table.

struct WenigerCache1F2{T}
    C::Vector{T}
    Cρ::Vector{T}
    C1::Vector{T}
    C2::Vector{T}
    C3::Vector{T}
    P::Vector{T}
    Q::Vector{T}
    N::Vector{T}
    ΔN::Vector{T}
    ΔNold::Vector{T}
    D::Vector{T}
    ΔD::Vector{T}
    ΔDold::Vector{T}
    R::Vector{T}
end

function WenigerCache1F2(T)
    p, q = 1, 2
    r = 5
    ρ = 3

    return WenigerCache1F2{T}(
        zeros(T, r), zeros(T, ρ+2), zeros(T, ρ+1), zeros(T, ρ+2), zeros(T, ρ+2), zeros(T, ρ+2),
        zeros(T, q+2), zeros(T, r+1), zeros(T, r+1), zeros(T, r+1), 
        zeros(T, r+1), zeros(T, r+1), zeros(T, r+1), zeros(T, r+1))
end

function Base.empty!(cache::WenigerCache1F2)
    fill!(cache.C, 0)
    fill!(cache.Cρ, 0)
    fill!(cache.C1, 0)
    fill!(cache.C2, 0)
    fill!(cache.C3, 0)
    fill!(cache.P, 0)
    fill!(cache.Q, 0)
    fill!(cache.N, 0)
    fill!(cache.ΔN, 0)
    fill!(cache.ΔNold, 0)
    fill!(cache.D, 0)
    fill!(cache.ΔD, 0)
    fill!(cache.ΔDold, 0)
    fill!(cache.R, 0)
end

function weniger1F2(α::T, β::AbstractVector{T}, z::T, cache; kmax::Int = 100_000) where T
    absα = abs(α)
    if norm(z) < eps(real(T)) || norm(prod(α)) < eps(prod(absα))
        return one(T)
    end
    γ = T(3)/2
    ζ = inv(z)
    p = 1
    q = 2

    empty!(cache)

    C,Cρ,C1,C2,C3,P,Q,N,ΔN,ΔNold,D,ΔD,ΔDold,R = (cache.C, cache.Cρ, cache.C1, cache.C2, cache.C3, 
        cache.P, cache.Q, cache.N, cache.ΔN, cache.ΔNold, cache.D, cache.ΔD, cache.ΔDold, cache.R)

    r = 5
    ρ = 3
    C[1] = one(T)
    Cρ[ρ+2] = one(T)
    for s in ρ:-1:0
        Cρ[s+1] = -(s+1)*Cρ[s+2]/(ρ+1-s)
    end
    C2[ρ+2] = one(T)/pochhammer(γ-ρ-2, ρ+2)
    C3[ρ+1] = one(T)/pochhammer(γ-ρ-1, ρ+2)
    C3[ρ+2] = one(T)/pochhammer(γ-ρ, ρ+2)
    t = γ
    for j in 1:p
        t *= α[j]+1
    end
    P[1] = t
    err = abs(γ)
    for j in 1:p
        err *= absα[j]+1
    end
    t = one(T)
    for j in 1:q
        t *= β[j]+1
    end
    Q[1] = 2t
    N[r+1] = prod(β)*ζ/prod(α)/(γ-1)
    ΔN[r] = N[r+1]/pochhammer(γ-ρ-1, ρ)
    D[r+1] = prod(β)*ζ/prod(α)/(γ-1)
    ΔD[r] = D[r+1]/pochhammer(γ-ρ-1, ρ)
    R[r+1] = N[r+1]/D[r+1]
    k = 0
    @inbounds while k < r || (k < kmax && errcheck(R[r], R[r+1], 10eps(real(T))))
        for j in 1:r
            N[j] = N[j+1]
            D[j] = D[j+1]
            R[j] = R[j+1]
        end
        t1 = zero(T)
        for j in 0:min(k, q+1)
            t1 += C[j+1]*Q[j+1]*ΔN[r-j]
        end
        if k ≤ ρ
            for j in 0:k
                t2 = one(T)
                for i in 1:q
                    t2 *= β[i]+j+1
                end
                t2 *= j+2
                t1 += C[j+1]*(-one(T))^(k-j)*pochhammer(j+γ, k-ρ-1)*t2
            end
        end
        t2 = zero(T)
        s2 = zero(T)
        for s in max(0, ρ+1-k):ρ
            s2 += Cρ[s+1]*C1[s+1]*(N[r-ρ+s]+(γ+k-ρ+s-2)*N[r-ρ+s-1])
        end
        s2 += (γ+k-1)*N[r]/pochhammer(γ+2k-ρ-1, ρ+2)
        t2 += P[1]*s2
        s2 = zero(T)
        for s in max(0, ρ+1-k):ρ+1
            s2 += Cρ[s+1]*C2[s+1]*(γ+2k-2ρ+2s-3)*N[r-ρ+s-1]
        end
        ΔNold[r+1] = s2
        for j in 1:min(k, p+1)
            t2 += C[j+1]*P[j+1]*(ΔNold[r+2-j] + ΔNold[r+1-j])
        end
        N[r+1] = ζ*t1-t2
        t1 = zero(T)
        for j in 0:min(k, q+1)
            t1 += C[j+1]*Q[j+1]*ΔD[r-j]
        end
        t2 = zero(T)
        s2 = zero(T)
        for s in max(0,ρ+1-k):ρ
            s2_c1 = (D[r-ρ+s]+(γ+k-ρ+s-2)*D[r-ρ+s-1])
            s2 += Cρ[s+1]*C1[s+1] * s2_c1
        end
        s2 += (γ+k-1)*D[r]/pochhammer(γ+2k-ρ-1, ρ+2)
        t2 += P[1]*s2
        s2 = zero(T)
        for s in max(0,ρ+1-k):ρ+1
            s2 += Cρ[s+1]*C2[s+1]*(γ+2k-2ρ+2s-3)*D[r-ρ+s-1]
        end
        ΔDold[r+1] = s2
        for j in 1:min(k, p+1)
            t2 += C[j+1]*P[j+1]*(ΔDold[r+2-j] + ΔDold[r+1-j])
        end
        D[r+1] = ζ*t1-t2
        if norm(P[1]) < eps(err)
            return N[r+1]/D[r+1]
        end
        N[r+1] /= P[1]/pochhammer(γ+2k-ρ-1, ρ+2)
        D[r+1] /= P[1]/pochhammer(γ+2k-ρ-1, ρ+2)
        R[r+1] = N[r+1]/D[r+1]
        s1 = zero(T)
        for s in max(0, ρ-k):ρ+1
            s1 += Cρ[s+1]*C3[s+1]*(γ+2k-2ρ+2s-1)*N[r-ρ+s]
        end
        ΔN[r+1] = s1
        s1 = zero(T)
        for s in max(0, ρ-k):ρ+1
            s1 += Cρ[s+1]*C3[s+1]*(γ+2k-2ρ+2s-1)*D[r-ρ+s]
        end
        ΔD[r+1] = s1
        k += 1
        for j in min(k, p+1):-1:0
            ΔNold[r-j] = (γ+2k-ρ-j-4)*ΔNold[r-j+1] + (k-j-1)*ΔNold[r-j]
        end
        for j in min(k, q+1):-1:0
            ΔN[r-j] = (γ+2k-ρ-j-2)*ΔN[r-j+1] + (k-j)*ΔN[r-j]
        end
        for j in min(k, p+1):-1:0
            ΔDold[r-j] = (γ+2k-ρ-j-4)*ΔDold[r-j+1] + (k-j-1)*ΔDold[r-j]
        end
        for j in min(k, q+1):-1:0
            ΔD[r-j] = (γ+2k-ρ-j-2)*ΔD[r-j+1] + (k-j)*ΔD[r-j]
        end
        for j in min(k, ρ):-1:1
            C[j+1] += C[j]
        end
        if k ≤ ρ+1
            for s in max(0, ρ+1-k):ρ
                C1[s+1] = pochhammer(k-ρ+s, ρ+1-s)/pochhammer(γ+2k-2ρ+s-2, ρ+2)
            end
            for s in max(0, ρ+1-k):ρ+1
                C2[s+1] = pochhammer(k-ρ+s, ρ+1-s)/pochhammer(γ+2k-2ρ+s-3, ρ+2)
            end
        else
            for s in 0:ρ
                C1[s+1] *= k/(k-one(T)-ρ+s)*(γ+2k-2ρ+s-4)/(γ+2k-ρ+s-2)*(γ+2k-2ρ+s-3)/(γ+2k-ρ+s-1)
            end
            for s in 0:ρ+1
                C2[s+1] *= k/(k-one(T)-ρ+s)*(γ+2k-2ρ+s-5)/(γ+2k-ρ+s-3)*(γ+2k-2ρ+s-4)/(γ+2k-ρ+s-2)
            end
        end
        if k ≤ ρ
            for s in max(0, ρ-k):ρ+1
                C3[s+1] = pochhammer(k-ρ+s+1, ρ+1-s)/pochhammer(γ+2k-2ρ+s-1, ρ+2)
            end
        else
            for s in max(0, ρ-k):ρ+1
                C3[s+1] *= (k+one(T))/(k-ρ+s)*(γ+2k-2ρ+s-3)/(γ+2k-ρ+s-1)*(γ+2k-2ρ+s-2)/(γ+2k-ρ+s)
            end
        end
        t = γ+k
        for j in 1:p
            t *= α[j]+k+1
        end
        err = abs(γ)+k
        for j in 1:p
            err *= absα[j]+k+1
        end
        for j in 2:p+2
            s = t - P[j-1]
            P[j-1] = t
            t = s
        end
        P[p+2] = t
        t = one(T)
        for j in 1:q
            t *= β[j]+k+1
        end
        t *= k+2
        for j in 2:q+2
            s = t - Q[j-1]
            Q[j-1] = t
            t = s
        end
        Q[q+2] = t
    end
    return isfinite(R[r+1]) ? R[r+1] : R[r]
end
