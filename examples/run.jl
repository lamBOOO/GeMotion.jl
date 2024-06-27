using GFluxx
using LineSearches: BackTracking, StrongWolfe


# Define boundary conditions
wave(x) = sin(pi*x[1])
basak_uniform = (;T_diri_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"], T_diri_expressions=[1.0,0.0,0.0,0.5,0.5])
basak_wave = (;T_diri_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"], T_diri_expressions=[wave,0.0,0.0,0.0,0.0])
turan = (;T_diri_tags=["leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint", "toprightpoint"], T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0])

# out = GFluxx.simulate(Pr=0.7, Ra=1E3, n=1.0, 50, linesearch=BackTracking(), levels=(;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))); turan...)

n_elems1 = 20
cases=
[
  # # turan
  # [1E3, 1E4, 0.6, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([i for i=1:2:13] |> x->vcat(x,-x))), turan],
  # [1E3, 1E4, 1.0, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:4.5],[5.0]) |> x->vcat(x,-x))), turan],
  # [1E3, 1E4, 1.8, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([i for i=0.1:0.2:1.1] |> x->vcat(x,-x))), turan],
  # [1E3, 1E5, 0.6, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([4,12,20,26,28,29] |> x->vcat(x,-x))), turan],
  # [1E3, 1E5, 1.0, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=1:2:11],[10]) |> x->vcat(x,-x))), turan],
  # [1E3, 1E5, 1.8, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([i for i=0.2:0.4:3.0] |> x->vcat(x,-x))), turan],
  # [1E3, 1E6, 0.6, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([10,30,50,60,65] |> x->vcat(x,-x))), turan],
  # [1E3, 1E6, 1.0, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=2:4:18],[19]) |> x->vcat(x,-x))), turan],
  # [1E3, 1E6, 1.8, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:5.5],[6]) |> x->vcat(x,-x))), turan],
]

outs = []
for case in cases
  out = GFluxx.simulate(Pr=case[1], Ra=case[2], n=case[3], n_elems=case[4], linesearch=case[5], levels=case[6]; case[7]...)
  push!(outs, out)
end
