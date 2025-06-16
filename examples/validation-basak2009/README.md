# Validation

- Fig 7-8 of:
```
Heat flow analysis for natural convection within trapezoidal enclosures based on heatline concept
Tanmay Basaka, S. Roy b, I. Pop
```
- The cases where temperature BCs jump can not be reproduced since Basak uses a Dirichlet value there, that is based on the average Nusselt number. However, the average Nusselt number uses the temperaure gradient, which goes to infinity and doesn't seem to convergence. THerefore, this approach yields different reults that are not simply off by a shifting or scaling but also deviate in the shape of the curves. A validation with Paraview showed that the aproach from Basak yields different heatlines when using Paraview to plot "grad(T) + U * T". Therefore, we stick with our approach.
