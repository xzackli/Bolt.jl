# Bolt

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xzackli.github.io/Bolt.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/Bolt.jl/dev)
[![Build Status](https://github.com/xzackli/Bolt.jl/workflows/CI/badge.svg)](https://github.com/xzackli/Bolt.jl/actions)
[![Coverage](https://codecov.io/gh/xzackli/Bolt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xzackli/Bolt.jl)

⚡ Bolt is a pure-Julia integrator for the Boltzmann equations in cosmology. In particular, it can accurately compute the gradient of the CMB power spectrum with respect to parameters using forward-mode automatic differentiation.

```julia
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
```

```
Result Comparison: 191.09327614330596 191.1227973108751
Simple Finite Difference:
  993.477 ms (1470317 allocations: 88.50 MiB)
ForwardDiff:
  1.333 s (2047395 allocations: 144.63 MiB)
191.1227973108751
```