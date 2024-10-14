<p align="center">
  <img height="100" src="media/logo_text.png" alt="Resume application project app icon">
<p>

<div align="center">
<img src="https://github.com/lamBOOO/GenMatFlow.jl/actions/workflows/test.yml/badge.svg" alt="Testing workflow badge"/>
</div>

# GenMatFlow.jl

- Navier-Stokes-Fourier Solver for generalized material laws, including:
  - Newtonian fluids: $\boldsymbol{\sigma} = \boldsymbol{D}(\boldsymbol{u}) = \frac12 (\boldsymbol{\nabla} \boldsymbol{u} + \boldsymbol{\nabla} \boldsymbol{u}^T)$
  - Non-Newtonian fluids with power-law: $\boldsymbol{\sigma} = K {\left(2 \boldsymbol{D} \boldsymbol{\colon}  \boldsymbol{D} \right)}^{\frac{1-n}{2}} \boldsymbol{D}(\boldsymbol{u})$
- Discretization using Finite Elements in [Gridap.jl](https://github.com/gridap/Gridap.jl)

## Installation

- Clone the repository and open the folder
```bash
git clone git@github.com:lamBOOO/GenMatFlow.jl.git
cd GenMatFlow.jl
```

- Install all Julia dependencies
```bash
julia --project -e 'import Pkg; Pkg.instantiate()'
```

- Run examples by navigating to the folder and execute the examples from the shell:
```bash
cd TODO
```
