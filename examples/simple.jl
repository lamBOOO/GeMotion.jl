using GeMotion
using LineSearches: BackTracking
using Gridap

# Setup model
# (P3)-(L6)-(P4)
#  |         |
# (L7)      (L8)
#  |         |
# (P1)-(L5)-(P2)
model = CartesianDiscreteModel((0, 1, 0, 1), (50, 50))
labels = get_face_labeling(model)

out = GeMotion.simulate(
  name="simple",
  Pr=0.7,
  Ra=1E3,
  n=1.0,
  model=model,
  jac_scaling=1,
  nlsolver_opts=(;
    show_trace=true,
    method=:newton,
    linesearch=BackTracking(),
    ftol=1E-8,
    xtol=1E-10
  );
  nlsolver_custom_init_guess=[],
  nlsolver_init_guess_type=:zero,
  (;
    T_diri_tags=[5, 7, 8, 1, 2],
    T_diri_expressions=[1.0, 0.0, 0.0, 0.5, 0.5],
    T_natural_tags=[3, 4, 6],
    V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
  )...
)

GeMotion.plot_all_unitsquare(
  out.psih, out.Pih, out.Th, out.uh, model, "simple",
  (;
    T=[0.1 * i for i = 1:10],
    psi=([0.01, 0.05, 0.1, 0.15, 0.18] |> x -> vcat(x, -x)),
    Pi=100,
    Sth=5,
    Sfl=5
  );
) |> display
