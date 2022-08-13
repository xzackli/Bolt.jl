using Bolt, Test
using DoubleFloats

const TOL = 1e-13  # tolerance for these moments
# if a higher tolerance is desired, the fault is in Bessels.jl's approximations.
# SpecialFunctions offers a more precise asymptotic behavior (1-2 digits at the end)

# for moderate-sizes arguments, both Weniger and asymptotic Lommel can work
@testset "moments: Bessel J, moderate argument" begin
    refs = (
        big"1.02810300641345785268082956068118",
        big"8.16960159665387499951281402709868",
        big"1148.91466190490173313233835015677")

    for (α, ref) ∈ zip( (0, 1, 2), refs) 
        T = Float64
        x, ν = T(200.0), T(2.5)
        α_minus_half = T(α - 1//2)  # Lommel-based stuff takes α-1/2 as input
        I1 = Bolt.J_moment_asymp(x, ν, α_minus_half)
        prefactor = Bolt.Bolt.J_moment_asymp_prefactor(T, ν, α)
        I2 = Bolt.J_moment_asymp_nu_five_halves(x, α_minus_half, prefactor)

        # Weniger does not perform well with very large arguments, check convergence
        T = Double64  # use more precision
        x, ν = T(x), T(ν)
        cache = Bolt.WenigerCache1F2(T)
        I3 = Bolt.J_moment_weniger_₁F₂(x, ν, α, cache)

        @test abs(1 - I1 / ref) < TOL
        @test abs(1 - I2 / ref) < TOL
        @test abs(1 - I3 / ref) < TOL
    end
end

## for small arguments (x < 50), we use Weniger + high precision floats
@testset "moments: Bessel J, small argument" begin
    T = Double64
    x, ν, α = T(10.0), T(2.5), T(1.5)
    cache = Bolt.WenigerCache1F2(T)
    I3 = Bolt.J_moment_weniger_₁F₂(x, ν, α, cache)

    ref = big"-0.98904817846028826228408967797229"
    @test abs(1 - I3 / ref) < 1e-15

    x, ν, α = T(0.1), T(2.5), T(-1//2)
    I3 = Bolt.J_moment_weniger_₁F₂(x, ν, α, cache)
    ref = big"0.00001772317062480308"
    @test abs(1 - I3 / ref) < TOL
end

##
@testset "moments: spherical Bessel J, large argument quadrupole" begin
    refs = (
        big"0.78787782588093298580386838516258",
        big"2.50028713446521582908074317765179",
        big"105.635871208275569897958010677779")
        

    for (k, ref) ∈ zip( (0, 1, 2), refs) 
        T = Float64
        x, ν = T(200.0), T(2)
        I1 = Bolt.sph_j_moment_asymp(x, ν, k)
        prefactor = Bolt.sph_j_moment_asymp_prefactor(T, ν, k)
        I2 = Bolt.sph_j_moment_asymp_nu_2(x, k, prefactor)

        # # Weniger does not perform well with very large arguments, check convergence
        T = Double64  # use more precision
        x, ν = T(x), T(ν)
        cache = Bolt.WenigerCache1F2(T)
        I3 = Bolt.sph_j_moment_weniger_₁F₂(x, ν, k, cache)

        @test abs(1 - I1 / ref) < TOL
        @test abs(1 - I2 / ref) < TOL
        @test abs(1 - I3 / ref) < TOL
    end
end

##
@testset "moments: spherical Bessel J, small argument quadrupole" begin
    T = Double64
    x, ν, α = T(0.1), T(2), T(1)
    cache = Bolt.WenigerCache1F2(T)
    I3 = Bolt.sph_j_moment_weniger_₁F₂(x, ν, α, cache)

    ref = big"1.66587318119689113603548520948514e-6"
    @test abs(1 - I3 / ref) < TOL

    x, ν, α = T(0.1), T(2), T(0)
    I3 = Bolt.sph_j_moment_weniger_₁F₂(x, ν, α, cache)
    ref = big"0.0000222127003021204915961147418458393"
    @test abs(1 - I3 / ref) < TOL

    x, ν, α = T(0.01), T(2), T(2)
    I3 = Bolt.sph_j_moment_weniger_₁F₂(x, ν, α, cache)
    ref = big"1.33332653062694211665894244946280e-12"
    @test abs(1 - I3 / ref) < TOL
end

## j₃ moments
@testset "moments: spherical Bessel J, large argument octupole" begin
    refs = (big"0.662361624450783328243409501185093",
        big"1.49770949285120051392753792923852",
        big"-163.183648420458825380490472539234")

    for (k, ref) ∈ zip( (0, 1, 2), refs) 
        T = Float64
        x, ν = T(200.0), T(3)
        I1 = Bolt.sph_j_moment_asymp(x, ν, k)
        prefactor = Bolt.sph_j_moment_asymp_prefactor(T, ν, k)
        I2 = Bolt.sph_j_moment_asymp_nu_3(x, k, prefactor)

        # # Weniger does not perform well with very large arguments, check convergence
        T = Double64  # use more precision
        x, ν = T(x), T(ν)
        cache = Bolt.WenigerCache1F2(T)
        I3 = Bolt.sph_j_moment_weniger_₁F₂(x, ν, k, cache)

        @test abs(1 - I1 / ref) < TOL
        @test abs(1 - I2 / ref) < TOL
        @test abs(1 - I3 / ref) < TOL
    end
end

## interpolation tests
@testset "interpolator: spherical Bessel J, large argument quadrupole" begin
    refs = (
        big"0.78787782588093298580386838516258",
        big"2.50028713446521582908074317765179",
        big"105.635871208275569897958010677779")

    order = 3
    ν = 2
    itp = Bolt.sph_bessel_interpolator(ν, order, 0.0, 1.6e4, 2_000_000)
    
    moms = itp(200.0)
    for i in 1:3
        @test abs(1 - moms[i] / refs[i]) < 1e-12  # interpolation error hurts!
    end
end

@testset "interpolator: spherical Bessel J, large argument octupole" begin
    refs = (big"0.662361624450783328243409501185093",
        big"1.49770949285120051392753792923852",
        big"-163.183648420458825380490472539234")

    order = 3
    ν = 3
    itp = Bolt.sph_bessel_interpolator(ν, order, 0.0, 1.6e4, 2_000_000)
    
    moms = itp(200.0)
    for i in 1:3
        @test abs(1 - moms[i] / refs[i]) < 1e-12  # interpolation error hurts!
    end
end

@testset "interpolator: spherical Bessel J, small argument octupole" begin
    refs = (2.380070697034449e-7,1.904006180460421e-8,1.586640331877486e-9)

    order = 3
    ν = 3
    itp = Bolt.sph_bessel_interpolator(ν, order, 0.0, 1.6e4, 2_000_000)
    
    moms = itp(0.1)
    for i in 1:3
        @test abs(moms[i] - refs[i]) < 1e-12  # interpolation error hurts!
    end
end

@testset "interpolator: spherical Bessel J, 4th order, moderate arg" begin
    refs = (big"0.00229425677577922706134041736035369", big"0.00183049831049510967778943596015052",
        big"0.00152235376642692867152374927665023", big"0.00130283667913981697274149102492971")

    order = 4
    ν = 3
    itp = Bolt.sph_bessel_interpolator(ν, order, 0.0, 1.6e4, 2_000_000)
    moms = itp(1.0)
    for i in 1:4
        @test abs(moms[i] - refs[i]) < 1e-12  # interpolation error hurts!
    end
end

@testset "interpolator: spherical Bessel J, 4th order, asymptotic arg" begin
    refs = (big"0.667496353968182420081865144062472", big"3.18644086496078566989171654620754",
        big"838.803790872929501155031335384019",big"831383.108409725479693899756952862")
    order = 4
    ν = 3
    itp = Bolt.sph_bessel_interpolator(ν, order, 0.0, 1.6e4, 2_000_000)
    moms = itp(1000.0)
    for i in 1:4
        @test abs(1 - moms[i] / refs[i]) < 1e-12  # interpolation error hurts!
    end
end

##
@testset "3rd order filon integration using the interpolator" begin

    ν, order = 3, 3
    itp = Bolt.sph_bessel_interpolator(ν, order, 0.0, 1.6e4, 2_000_000)
    
    f(x) = 3x^2 - 0.2x + 4
    
    @test abs(Bolt.integrate_sph_bessel_filon(4., -0.2, 6.0, 10.0, 0., 2., itp) - 0.365287615501162668736682652658444) < TOL
    @test abs(Bolt.integrate_sph_bessel_filon(4., -0.2, 6.0, 10.0, 0., 1., itp) - 0.219009396999160523658045931736310) < TOL
    @test abs(Bolt.integrate_sph_bessel_filon(6.8, 5.8, 6.0, 10.0, 1., 2., itp) - 0.146278218502002145078636720922135) < TOL
    
    itp_ka = itp(10 * 1.)  # k*a
    s, _ = Bolt._loop_integrate_sph_bessel_filon(6.8, 5.8, 6.0, 10.0, 1., 2., itp, itp_ka)
    @test abs(s - 0.146278218502002145078636720922135) < TOL
end
