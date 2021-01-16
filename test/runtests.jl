using Bolt
using Test

@testset "Bolt.jl" begin
    # Write your tests here.

    # early_time_Xₑ(z2x(1587.4))
end

## hand-written ℋ derivative for testing
# function ℋ′(x, par::AbstractCosmoParams)
#     a = x2a(x)
#     return -H₀(par) * (2par.Ω_r + (par.Ω_b + par.Ω_m) * a - 2Ω_Λ(par) * a^4) /
#         (2 * a * √(par.Ω_r + (par.Ω_b + par.Ω_m) * a + Ω_Λ(par) * a^4))
# end
