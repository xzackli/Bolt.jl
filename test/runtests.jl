using Bolt
using Test
using DelimitedFiles
using LinearAlgebra
using ForwardDiff

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
    𝕡 = CosmoParams(Σm_ν=0.0, N_ν=3.0, Ω_r=5.042e-5)
    bg = Background(𝕡)
    𝕣 = Bolt.RECFAST(bg=bg, OmegaB=𝕡.Ω_b, Yp=𝕡.Y_p, OmegaG=𝕡.Ω_r, Tnow=2.725)
    xe_bespoke, Tmat_bespoke = Bolt.recfast_xe(𝕣; Nz=1000, zinitial=10000., zfinal=0.)
    #change to only test pre-reion (z≧50)
    @test all(abs.(Xe_fort[1:end-5] .- xe_bespoke[1:end-5]) .< 1e-5)
end

##

#Diff tests for bg and ih+𝕣 #FIXME these can probably just be one test?
#bg
@testset "bg_fwddiff" begin
    function fbg(Ω_b::DT) where DT
       𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
       bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
       return bg.η(-5)
    end
    fbg(0.046)
    Δ = 1e-3
    (fbg(0.046+ Δ) - fbg(0.046 - Δ)) / 2Δ
    @test (((fbg(0.046+ Δ) - fbg(0.046 - Δ)) / 2Δ - ForwardDiff.derivative(fbg, 0.046)) .< 1e-5)
end

# ih with recfast
@testset "ih_fwddiff" begin
    function fih(Ω_b::DT) where DT
       𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
       bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
       𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
       #TODO?: Need to supply all three relevant cosmo params to recfast to avoid dual problem
       ih = IonizationHistory(𝕣, 𝕡, bg)
       return ih.csb²(0.)
    end
    fih(0.046)
    Δ = 1e-3
    (fih(0.046+ Δ) - fih(0.046 - Δ)) / 2Δ
    @test (((fih(0.046+ Δ) - fih(0.046 - Δ)) / 2Δ - ForwardDiff.derivative(fih, 0.046)) .< 1e-5)
end
##
