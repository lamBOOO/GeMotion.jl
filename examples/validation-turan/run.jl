using GeMotion
using Gridap
using LineSearches: BackTracking, StrongWolfe


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
model_small = CartesianDiscreteModel((0, 1, 0, 1), (20, 20))  # see later
labels_small = get_face_labeling(model_small)
add_tag_from_tags!(labels_small, "botleftpoint", [1,])
add_tag_from_tags!(labels_small, "botrightpoint", [2,])
add_tag_from_tags!(labels_small, "topleftpoint", [3,])
add_tag_from_tags!(labels_small, "toprightpoint", [4,])
add_tag_from_tags!(labels_small, "botline", [5,])
add_tag_from_tags!(labels_small, "topline", [6,])
add_tag_from_tags!(labels_small, "leftline", [7,])
add_tag_from_tags!(labels_small, "rightline", [8,])
add_tag_from_tags!(labels_small, "all", [1, 2, 3, 4, 5, 6, 7, 8])

# Define boundary conditions
turan = (;
  T_diri_tags = [
    "leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint",
    "toprightpoint"
    ],
  T_diri_expressions = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
)

nlsolver_opts = (;
  show_trace = true,
  method = :newton,
  ftol = 1E-8,
  xtol = 1E-50,
  linesearch = BackTracking(),
)

cases =
[
  # turan
  [1E3, 1E4, 0.6, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([i for i=1:2:13] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E4, 1.0, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:4.5],[5.0]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E4, 1.8, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([i for i=0.1:0.2:1.1] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E5, 0.6, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([4,12,20,26,28,29] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E5, 1.0, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=1:2:11],[10]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E5, 1.8, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([i for i=0.2:0.4:3.0] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E6, 0.6, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([10,30,50,60,65] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E6, 1.0, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=2:4:18],[19]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E6, 1.8, model, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:5.5],[6]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
]

# Make faster when in GitHub Actions environment
if haskey(ENV, "GITHUB_ACTIONS")
# if true
  cases = cases[1:1]
  cases[1][4] = model_small
  @info "Changed values for CI:"
  @show cases
end

outs = []
for (i,case) in enumerate(cases)
  out = GeMotion.simulate(name="$i", Pr=case[1], Ra=case[2], n=case[3], model=case[4], nlsolver_opts=case[5]; case[7]...)
  push!(outs, out)
  out2 = GeMotion.plot_all_unitsquare(
    out.psih, out.Th, out.uh, case[4], "$i", case[6]
  )
end
