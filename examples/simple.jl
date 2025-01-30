using GeMotion
using LineSearches: BackTracking
using Gridap

# Setup model
# labelling:
# 1-2-3-4 = botleftpoint-botrightpoint-topleftpoint-toprightpoint
# 5-6-7-8 = botline-topline-leftline-rightline
model = CartesianDiscreteModel((0, 1, 0, 1), (50, 50))
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "botleftpoint", [1,])
add_tag_from_tags!(labels, "botrightpoint", [2,])
add_tag_from_tags!(labels, "topleftpoint", [3,])
add_tag_from_tags!(labels, "toprightpoint", [4,])
add_tag_from_tags!(labels, "botline", [5,])
add_tag_from_tags!(labels, "topline", [6,])
add_tag_from_tags!(labels, "leftline", [7,])
add_tag_from_tags!(labels, "rightline", [8,])
add_tag_from_tags!(labels, "all", [1, 2, 3, 4, 5, 6, 7, 8])

out = GeMotion.simulate(
  name="simple",
  Pr=0.7,
  Ra=1E3,
  n=1.0,
  model=model,
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

GeMotion.plot_all_unitsquare(
  out.psih, out.Th, out.uh, model, "simple",
  (;
    T=[0.1 * i for i = 1:10],
    psi=([0.01, 0.05, 0.1, 0.15] |> x -> vcat(x, -x)),
    Sth=5,
    Sfl=5
  );
)
