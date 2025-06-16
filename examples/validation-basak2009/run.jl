using GeMotion
using LineSearches: BackTracking, StrongWolfe, Static, HagerZhang, MoreThuente
using Makie
using Gridap
using Gridap.Arrays  # for "return_cache"
using CSV, DataFrames

# Setup model
# (P3)-(L6)-(P4)
#  |         |
# (L7)      (L8)
#  |         |
# (P1)-(L5)-(P2)
model = CartesianDiscreteModel(
  (0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (50, 50) : (100, 100)
)

# solver settings
nlsolver_opts = (;
  show_trace=true,
  method=:newton,
  linesearch=BackTracking(),
  ftol=1E-8,
  xtol=1E-50,
  iterations=200,
)

# BCs
uniform = (;
  T_diri_tags=[5, 7, 8, 1, 2],
  T_diri_expressions=[1.0, 0.0, 0.0, 0.5, 0.5],
  T_natural_tags=[3, 4, 6],
  V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
)
wave = (;
  T_diri_tags=[5, 7, 8, 1, 2],
  T_diri_expressions=[x -> sin(pi * x[1]), 0.0, 0.0, 0.0, 0.0],
  T_natural_tags=[3, 4, 6],
  V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
)
nonuniform_side = (;
  T_diri_tags=[5, 7, 8, 1, 2],
  T_diri_expressions=[1.0, x -> 1 - x[2], x -> 1 - x[2], 1.0, 1.0],
  T_natural_tags=[3, 4, 6],
  V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
)

common_lvls = (; T=[0.1 * i for i = 1:10])

paramlist = [
# Pr, Ra, n, psi_lvls, Pi_lvls, Sth_lvls, Sfl_lvls
(0.026, 1E3, 1.0, [0.001; 0.01; 0.05:0.05:0.15; 0.18],
  [0.01; 0.05; 0.1; 0.15; 0.25; 0.35; 0.55; 0.7], 5, 5, uniform,
  "basak2009_fig3_diff"),
(0.026, 1E5, 1.0, [0.01; 0.05; 0.1; 0.5; 1.5; 3.5; 5.5; 6],
  [0.05; 0.4; 0.6; 0.7; 0.8; 1; 1.2], 5, 5, uniform,
  "basak2009_fig4_diff"),
(0.7, 1E5, 1.0, [0.1; 0.8; 2; 4.0; 7; 10.0; 13],
  [0.1; 0.5; 1.5; 2.5; 3; 3.5; 4; 5.5; 6.0], 5, 5, uniform,
  "basak2009_fig5_diff"),
(1000, 1E5, 1.0, [0.1; 1; 4; 7; 10; 13; 15],
  [0.5; 1.5; 2.5; 3.5; 4.5; 5.5; 6.5; 7], 5, 5, uniform,
  "basak2009_fig6_diff"),
(0.026, 1E5, 1.0, [0.05; 0.5; 2; 3.5; 5],
  [0.05; 0.4; 0.6; 0.7; 0.8; 1; 1.2], 5, 5, wave,
  "basak2009_fig7_same"),
(0.7, 1E5, 1.0, [0.1; 0.8; 2; 4.0; 7; 10.0; 12],
  [0.1; 0.5; 1.0; 1.5; 2.0; 2.5; 3; 4; 4.5], 5, 5, wave,
  "basak2009_fig8_same"),
(0.015, 1E3, 1.0, [0.001; 0.005; 0.01; 0.02; 0.04; 0.06],
  [0.005; 0.04; 0.1; 0.15; 0.25], 5, 5, nonuniform_side,
  "basak2011_fig2_diff"),
(0.015, 1E4, 1.0, [0.001; 0.005; 0.01; 0.04; 0.1; 0.3; 0.6; 0.8],
  [0.01; 0.05; 0.1; 0.2; 0.3; 0.5; 0.6], 5, 5, nonuniform_side,
  "basak2011_fig3_diff"),
(0.7, 1E5, 1.0, [0.1; 1; 2; 3; 5; 6],
  [0.1; 0.5:0.5:4.5], 5, 5, nonuniform_side,
  "basak2011_fig5_diff"),
]
function mkcase(Pr, Ra, n, psi_vec, Pi_vec, Sth_vec, Sfl_vec, bcs, name)
  [
    Pr, Ra, n, model, nlsolver_opts,
    (;
      common_lvls...,
      psi=typeof(psi_vec) == Int ? psi_vec : vcat(psi_vec, -psi_vec),
      Pi=typeof(Pi_vec) == Int ? Pi_vec : vcat(Pi_vec, -Pi_vec),
      Sth=Sth_vec, Sfl=Sfl_vec
    ), bcs, name
  ]
end
cases = [
  mkcase(Pr, Ra, n, psi_vec, Pi_vec, Sth_vec, Sfl_vec, bcs, name)
  for (Pr, Ra, n, psi_vec, Pi_vec, Sth_vec, Sfl_vec, bcs, name) in paramlist
]

# Run cases
outs = []
outs2 = []
for (i, case) in enumerate(cases)
  out2 = GeMotion.simulate(
    name="$(case[8])_$(case[1])_$(case[2])", Pr=case[1], Ra=case[2], n=case[3],
    model=case[4], jac_scaling=1, nlsolver_opts=case[5],
    nlsolver_custom_init_guess=[], nlsolver_init_guess_type=:zero, ; case[7]...
  )
  out = GeMotion.plot_all_unitsquare(
    out2.psih, out2.Pih, out2.Th, out2.uh, model, "$(case[8])_$(case[1])_$(case[2])", case[6]
  )
  push!(outs2, out2)
  push!(outs, (; Pr=case[1], Ra=case[2], out...))
end
