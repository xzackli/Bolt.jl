# Bolt

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/Bolt.jl/dev)
[![](https://img.shields.io/badge/source-github-blue)](https://github.com/xzackli/Bolt.jl) 

[![Build Status](https://github.com/xzackli/Bolt.jl/workflows/CI/badge.svg)](https://github.com/xzackli/Bolt.jl/actions)
[![codecov](https://codecov.io/gh/xzackli/Bolt.jl/branch/main/graph/badge.svg?token=NDj9hvOUkN)](https://codecov.io/gh/xzackli/Bolt.jl)

⚡ Bolt is a pure-Julia integrator for the Boltzmann equations in cosmology. It can accurately compute the gradient of the CMB power spectrum, with respect to cosmological parameters, using forward-mode automatic differentiation.

**Contributors**: Jamie Sullivan, Zack Li, Marius Millea

## Install

Bolt requires Julia 1.5+. To install, from the package prompt, run:

```
pkg> add https://github.com/xzackli/Bolt.jl
```

## Gallery

*A CMB temperature power spectrum and gradient from Bolt.jl.*
![](example_spectrum.png) ![](docs/src/example_spectrum.png)

*A linear matter power spectrum and gradient from Bolt.jl.*
![](example_linear_power_c.png) ![](docs/src/example_linear_power_c.png)


