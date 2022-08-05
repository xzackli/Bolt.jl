
using Test
using DoubleFloats

# for moderate-sizes arguments, both Weniger and asymptotic Lommel can work
@testset "moments: Bessel J, moderate argument" begin
    T = Float64
    x, ν, α = T(200.0), T(2.5), T(-0.5)
    
    α_minus_half = T(α - 1//2)  # Lommel-based stuff takes α-1/2 as input
    I1 = Bolt.J_moment_asymptotic_lommel(x, ν, α_minus_half)
    prefactor = Bolt.Bolt.J_moment_asymptotic_lommel_prefactor(T, ν, α)
    I2 = Bolt.J_moment_asymptotic_nu_five_halves(x, α_minus_half, prefactor)

    # Weniger does not perform well with very large arguments, check convergence
    T = Double64  # use more precision
    x, ν, α = T(x), T(ν), T(α)
    cache = Bolt.WenigerCache1F2(T)
    I3 = Bolt.J_moment_weniger_₁F₂(x, ν, α, cache)

    ref = big"0.62863555306932463883337671258395"
    @test abs(1 - I1 / ref) < 1e-15
    @test abs(1 - I2 / ref) < 1e-15
    @test abs(1 - I3 / ref) < 1e-15
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
    @test abs(1 - I3 / ref) < 1e-15
end
