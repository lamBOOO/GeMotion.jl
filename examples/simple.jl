using GeMotion
using LineSearches: BackTracking
using Gridap

out = GeMotion.simulate(
  name="simple",
  Pr=0.7,
  Ra=1E3,
  n=1.0,
  n_elems=50,
  nlsolver_opts=(;
    show_trace=true,
    method=:newton,
    linesearch=BackTracking(),
    ftol=1E-8,
    xtol=1E-10
  ),
  levels=(;
    T=[0.1 * i for i = 1:10],
    psi=([0.01, 0.05, 0.1, 0.15] |> x -> vcat(x, -x)),
    Sth=5,
    Sfl=5
  );
  (;
    T_diri_tags=[
      "leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint",
      "toprightpoint"
    ],
    T_diri_expressions=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
  )...
)

