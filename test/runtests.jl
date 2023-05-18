using Bolt
using Test
using DelimitedFiles
using LinearAlgebra
using ForwardDiff
using DataInterpolations
using Printf

include("testbessel.jl")

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
    z⃗, Xe_fort = recfastdata[:,1][1:end-5], recfastdata[:,2][1:end-5]
    𝕡 = CosmoParams(Σm_ν=0.0, N_ν=3.0, Ω_r=5.042e-5)
    bg = Background(𝕡)
    𝕣 = Bolt.RECFAST(bg; OmegaB=𝕡.Ω_b, Yp=𝕡.Y_p, OmegaG=𝕡.Ω_r, Tnow=2.725)
    rhist = Bolt.recfastsolve(𝕣)
    # xe_bespoke, Tmat_bespoke = Bolt.recfast_xe(𝕣; Nz=1000, zinitial=10000., zfinal=0.)
    #change to only test pre-reion (z≧50)
    @test all(abs.(Xe_fort .- Bolt.Xe.((rhist,), z⃗)) .< 1e-4)
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

## ih with recfast
@testset "ih_fwddiff" begin
    function fih(Ω_b::DT) where DT
       𝕡 = CosmoParams{DT}(Ω_b=Ω_b)
       bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=15)
       𝕣 = Bolt.RECFAST(bg; Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r, tol=1e-10)
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

@testset "nonu_class_comparison_1e-3" begin
    # bg/ion setup
    𝕡 = CosmoParams(Σm_ν=0.0)
    n_q=15
    logqmin,logqmax = -6,-1
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=n_q)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    ih = IonizationHistory(𝕣, 𝕡, bg)


    x_grid = bg.x_grid

    # Choose a k-mode to compare to saved class perturbations at
    k_options = ["p03", "p3", "1p0", #choose from k = [0.03h/Mpc, 0.3h/Mpc, 1.0h/Mpc]
                "p01", ] #newly computed k modes
    k_choice = k_options[1]
    #Read in CLASS perturbations
    #CLASS keys (for reference):
    #['k (h/Mpc)', 'd_g', 'd_b', 'd_cdm', 'd_ur', 'd_ncdm[0]', 'd_tot',
    #'phi', 'psi', 't_g', 't_b', 't_cdm', 't_ur', 't_ncdm[0]', 't_tot']
    retnf = open( @sprintf("data/zack_N_class_px_k%s_nofluid_nonu.dat",k_choice),"r" ) do datafile
    # an example that goes to early times -> retnf = open("./test/data/lowres_class_px_kp03_nofluid.dat","r") do datafile
        [parse.(Float64, split(line)) for line in eachline(datafile)]
    end
    #the second column is just a repeated k value, so remember it and delete col
    kclass = retnf[2][1] #read class k mode from file (in h/Mpc)
    # k = (bg.H₀*3e5/100)*kclass #get k in our units ->old value
    k = 𝕡.h * kclass  #get k in our units
    class_pxsnf = transpose(reduce(hcat,retnf[1:end .!= 2]))

    xhor = x_grid[argmin(abs.(k ./ (2π* bg.ℋ.(x_grid).*𝕡.h) .- 1))] #horizon crossing ish
    println("k = ", kclass," log10k = ", log10(kclass), " h/Mpc")

    #pert setup
    ℓᵧ=50
    ℓ_ν=50
    ℓ_mν=20
    reltol=1e-9
    abstol=1e-9
    pertlen = 2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+5
    results=zeros(pertlen,length(x_grid))
    ℳρ,ℳσ = zeros(length(x_grid)),zeros(length(x_grid)) #arrays for the massive neutrino integrated perts
    hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,n_q)
    #solve (with rsa)
    results_with_rsa = boltsolve_rsa(hierarchy; reltol=reltol, abstol=abstol)

    class_x = class_pxsnf[1,:][end:-1:1]

    itphibolt = CubicSpline((results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+1,:]), x_grid)
    itpphiclass = CubicSpline(class_pxsnf[7,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

    itdeltbbolt = CubicSpline((results_with_rsa[2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*n_q+4,:]), x_grid)
    itdeltbclass = CubicSpline(class_pxsnf[3,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

    itpgambolt = CubicSpline(-(results_with_rsa[1,:]*4)[1:end], x_grid)
    itpgamclass = CubicSpline(class_pxsnf[2,:][end:-1:1], class_pxsnf[1,:][end:-1:1])

    class_eta = bg.η.(class_x)


    TOL = 1e-3
    @test all(abs.(itphibolt.(class_x) ./ itpphiclass.(class_x) .- 1) .< TOL)
    @test all(abs.(-itdeltbbolt.(class_x) ./ itdeltbclass.(class_x) .- 1) .< TOL)
    @test all(abs.(itpgambolt.(class_x) ./ itpgambolt.(class_x) .- 1) .< TOL)
end

@testset "cell_camb_1e-1" begin
    # the usual
    𝕡 = CosmoParams()
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=15)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100) #FIXME be careful with this for derivatives...
    # source functions
    sf_t = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    sf_e = source_grid_P(𝕡, bg, ih, k_grid, BasicNewtonian())

    ℓs = 10:10:2500 #pretty coarse but ok - enough to see wiggles even at 20
    Cℓtt = cltt(ℓs, 𝕡, bg, ih, sf_t)
    Cℓte = clte(ℓs, 𝕡, bg, ih, sf_t,sf_e)
    Cℓee = clee(ℓs, 𝕡, bg, ih, sf_e)

    #camb results
    # class_Cℓs = readdlm("data/class_rough_ttteee_unlensed.dat")
    camb_Cℓs = readdlm("data/camb_rough_ttteee_unlensed.dat")

    #FIXME update
    #ratios
    # itptt = linear_interpolation(camb_Cℓs[1,:], camb_Cℓs[2,:])
    # itpte = linear_interpolation(camb_Cℓs[1,:], camb_Cℓs[3,:])
    # itpee = linear_interpolation(camb_Cℓs[1,:], camb_Cℓs[4,:])
    # #make same interpolation to ensure same interp error
    # itptt_b = linear_interpolation(ℓs, @.(ℓfac * Cℓtt))
    # itpte_b = linear_interpolation(ℓs, @.(ℓfac * Cℓte))
    # itpee_b = linear_interpolation(ℓs, @.(ℓfac * Cℓee))
    ℓfac = ℓs.*(ℓs.+1)
    itptt = LinearInterpolation(camb_Cℓs[2,:], camb_Cℓs[1,:])
    itpte = LinearInterpolation(camb_Cℓs[3,:], camb_Cℓs[1,:])
    itpee = LinearInterpolation(camb_Cℓs[4,:], camb_Cℓs[1,:])
    #make same interpolation to ensure same interp error
    itptt_b = LinearInterpolation(@.(ℓfac * Cℓtt), ℓs)
    itpte_b = LinearInterpolation(@.(ℓfac * Cℓte), ℓs)
    itpee_b = LinearInterpolation(@.(ℓfac * Cℓee), ℓs)

    #test line
    TOL = 1.1e-1
    #FIXME should really do this on the finer of the two grids but need to rerun...
    @test all(abs.(itptt_b.(ℓs)  ./ itptt.(ℓs) .- 1)  .< TOL)
    # all(abs.(itpte_b.(ℓs)  .- itpte.(ℓs) .- 1)  .< TOL) # zero-crossing will be an issue here, do diff check
    @test all(abs.(itpee_b.(ℓs)  ./ itpee.(ℓs) .- 1)  .< TOL)
end
