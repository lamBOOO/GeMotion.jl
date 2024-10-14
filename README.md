<p align="center">
  <img height="100" src="media/logo_textcentered.png" alt="Resume application project app icon">
<p>

<div align="center">
  <a href="https://github.com/lamBOOO/GeMotion.jl/actions"><img src="https://github.com/lamBOOO/GeMotion.jl/actions/workflows/test.yml/badge.svg" alt="Testing workflow badge"/></a>
  <a href="https://lambooo.github.io/GeMotion.jl/dev/"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev documentation badge"/></a>
  <a href="https://github.com/lamBOOO/GeMotion.jl/blob/master/LICENSE"><img src="https://img.shields.io/github/license/lamBOOO/GeMotion.jl.svg" alt="License badge"/></a>
</div>

# GeMotion.jl: A Navier-Stokes-Fourier Solver for generalized material laws

- Different material laws, including:
  - Newtonian fluids: $\boldsymbol{\sigma} = \boldsymbol{D}(\boldsymbol{u}) = \frac12 (\boldsymbol{\nabla} \boldsymbol{u} + \boldsymbol{\nabla} \boldsymbol{u}^T)$
  - Non-Newtonian fluids with power-law: $\boldsymbol{\sigma} = K {\left(2 \boldsymbol{D} \boldsymbol{\colon}  \boldsymbol{D} \right)}^{\frac{1-n}{2}} \boldsymbol{D}(\boldsymbol{u})$
- Discretization using Finite Elements in [Gridap.jl](https://github.com/gridap/Gridap.jl)
- Solve nonlinear systems using Newtons method


## Installation

- Clone the repository and open the folder
```bash
git clone git@github.com:lamBOOO/GeMotion.jl.git
cd GeMotion.jl
```

- Install all Julia dependencies
```bash
julia --project -e 'import Pkg; Pkg.instantiate()'
```

- Run examples by navigating to the folder and execute the examples from the shell:
```bash
cd examples/validation-basak
julia --project basak.jl
```
