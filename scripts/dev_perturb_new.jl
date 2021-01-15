using Bolt
using ForwardDiff

# z=0 metric perturbation
function Φ₀(Ω_b::DT) where DT
    par = CosmoParams{DT}(Ω_b=Ω_b)
    bg = Background(par)
    ih = IonizationHistory(SahaPeebles(), par, bg)
    hierarchy = Hierarchy(340bg.H₀, 8, par, bg, ih)
    metric_perturbation = boltsolve(hierarchy, BasicNewtonian(); reltol=1e-10)
    return metric_perturbation(0.0)[1]
end

Δ = 1e-6
print("Result Comparison: ",
    (Φ₀(0.046 + Δ) - Φ₀(0.046 - Δ)) / 2Δ, " ",
    ForwardDiff.derivative(Φ₀, 0.046), "\n")


using BenchmarkTools
print("Simple Finite Difference:\n")
@btime (Φ₀(0.046 + Δ) - Φ₀(0.046 - Δ)) / 2Δ
print("ForwardDiff:\n")
@btime ForwardDiff.derivative(Φ₀, 0.046)
