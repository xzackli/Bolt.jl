using Bolt
using Test
using DelimitedFiles
using LinearAlgebra

@testset "FFTLog" begin
    N = 64
    μ = 0
    q = 0.0
    r₀ = 1.0
    L = 8.0
    Nhalf = N ÷ 2
    n = range(-Nhalf,Nhalf,length=N)
    r = r₀ .* 10 .^ (n .* L ./ N )
    pl = Bolt.plan_fftlog(r, μ, q, 1.0; kropt=true)
    aₙ = r .^ (μ + 1) .* exp.(-r.^2 / 2)
    y = similar(r, ComplexF64)
    fftdata = readdlm("data/fftlog_example.txt", ' ', Float64, '\n')

    # test forward
    mul!(y, pl, aₙ)
    f_ref = fftdata[:,2]
    @test all(abs.(y .- f_ref) .< 1e-15)
    @test isapprox(y, f_ref)

    # test backward
    y2 = similar(r, ComplexF64)
    ldiv!(y2, pl, y)
    @test all(abs.(y2 .- aₙ) .< 1e-15)
end

##
@testset "RECFAST" begin
    recfastdata = readdlm("data/test_recfast_1.dat", ',', Float64, '\n', header=true)[1]
    z⃗, Xe_fort = recfastdata[:,1], recfastdata[:,2]
    𝕡 = CosmoParams(Σm_ν=0.0, N_ν=3.0)
    bg = Background(𝕡)
    𝕣 = Bolt.RECFAST(bg=bg, OmegaB=𝕡.Ω_b, Yp=𝕡.Y_p)
    xe_bespoke, Tmat_bespoke = Bolt.recfast_xe(𝕣; Nz=1000, zinitial=10000., zfinal=0.)
    @test all(abs.(Xe_fort .- xe_bespoke) .< 1e-5)
end
