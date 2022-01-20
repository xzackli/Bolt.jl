# Bolt

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xzackli.github.io/Bolt.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/Bolt.jl/dev)
[![Build Status](https://github.com/xzackli/Bolt.jl/workflows/CI/badge.svg)](https://github.com/xzackli/Bolt.jl/actions)
[![codecov](https://codecov.io/gh/xzackli/Bolt.jl/branch/main/graph/badge.svg?token=NDj9hvOUkN)](https://codecov.io/gh/xzackli/Bolt.jl)

⚡ Bolt is a pure-Julia integrator for the Boltzmann equations in cosmology. It can accurately compute the gradient of the CMB power spectrum, with respect to cosmological parameters, using forward-mode automatic differentiation.

It needs a fair bit more physics before it can be applied to modern data. There are some examples showing accurate reproduction of some figures in `examples/`.

I haven't spent any time optimizing Bolt for performance yet, and in particular the naively implemented source function integrations currently dominate the spectrum cost.

<p align="center">
<img width=80% src="docs/assets/example_spectrum.png">
</p>

*A CMB temperature power spectrum and gradient from Bolt.jl.*

<p align="center">
<img width=80% src="docs/assets/example_linear_power_c.png">
</p>

*A linear matter power spectrum and gradient from Bolt.jl.*

