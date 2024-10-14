using GeMotion
using LineSearches: BackTracking, StrongWolfe


# Define boundary conditions
turan = (;T_diri_tags=["leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint", "toprightpoint"], T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0])

nlsolver_opts = (;
  show_trace = true,
  method = :newton,
  ftol = 1E-8,
  xtol = 1E-50,
  linesearch = BackTracking(),
)

n_elems = 200
if haskey(ENV, "GITHUB_ACTIONS")
  # Lower the value when in GitHub Actions environment
  n_elems = 20  # New value for GitHub Actions
end

cases =
[
  # turan
  [1E3, 1E4, 0.6, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([i for i=1:2:13] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E4, 1.0, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:4.5],[5.0]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E4, 1.8, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([i for i=0.1:0.2:1.1] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E5, 0.6, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([4,12,20,26,28,29] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E5, 1.0, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=1:2:11],[10]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E5, 1.8, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([i for i=0.2:0.4:3.0] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E6, 0.6, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=([10,30,50,60,65] |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E6, 1.0, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=2:4:18],[19]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
  [1E3, 1E6, 1.8, n_elems, nlsolver_opts, (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:5.5],[6]) |> x->vcat(x,-x)), Sth=5, Sfl=10), turan],
]

outs = []
for (i,case) in enumerate(cases)
  out = GeMotion.simulate(name="$i", Pr=case[1], Ra=case[2], n=case[3], n_elems=case[4], nlsolver_opts=case[5], levels=case[6]; case[7]...)
  push!(outs, out)
end
