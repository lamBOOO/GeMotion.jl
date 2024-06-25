using GFluxx
using LineSearches: BackTracking, StrongWolfe


# Define boundary conditions
wave(x) = sin(pi*x[1])
basak_uniform = (;bcname="basak_uniform", T_diri_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"], T_diri_expressions=[1.0,0.0,0.0,0.5,0.5])
basak_wave = (;bcname="basak_wave", T_diri_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"], T_diri_expressions=[wave,0.0,0.0,0.0,0.0])
turan = (;bcname="turan", T_diri_tags=["leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint", "toprightpoint"], T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0])

# out = GFluxx.simulate(Pr=0.7, Ra=1E3, n=1.0, 50, linesearch=BackTracking(), levels=(;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))); turan...)

cases=
[
  [0.7, 1E3, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_uniform],
  [0.7, 5*1E3, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.15,0.5,1,1.3] |> x->vcat(x,-x))), basak_uniform],
  [0.7, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,13] |> x->vcat(x,-x))), basak_uniform],
  [0.1, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,4,7,9] |> x->vcat(x,-x))), basak_uniform],
  [1.0, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_uniform],
  [10.0, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_uniform],
  [0.015, 1E3, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_uniform],
  [0.7, 1E3, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_wave],
  [0.7, 5*1E3, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.15,0.5,1,1.3] |> x->vcat(x,-x))), basak_wave],
  [0.7, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,13] |> x->vcat(x,-x))), basak_wave],
  [0.1, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,4,7,9] |> x->vcat(x,-x))), basak_wave],
  [1.0, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_wave],
  [10.0, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_wave],
  [0.015, 1E3, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_wave],
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


# Nusselt number
begin
nn = 100
xs = 0 |> x->LinRange(0+x, 1-x, nn)
# update_theme!(
#   palette = (color=Makie.wong_colors(), marker = [:circle, :xcross, :pentagon, :diamond],),
#   ScatterLines = (cycle = Cycle([:color, :marker], covary = true),)
# )
update_theme!(
  palette = (color=Makie.wong_colors(), marker = [:circle, :xcross, :pentagon, :diamond],),
  Lines = (cycle = Cycle([:color, :marker], covary = true),)
)
f = Figure(
  # size = (500, 500)
  figure_padding = (0, 10, 0, 0)
)
ax = Axis(
  f[1, 1]
  ,aspect=300/180
  ,title="local Nusselt number Nu, bottom wall"
  ,xlabel="x"
  ,ylabel="y"
  ,xminorticksvisible = true
  ,yminorticksvisible = true
  ,xticks = LinearTicks(5)
  ,yticks = LinearTicks(3)
  ,limits = ((0.1, 0.9), (0,15))
)

map(outs[[1,3,6]]) do o
  cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
  lines!(xs, [evaluate!(cache, o.Nu, Gridap.Point([x,0])) for x in xs]
    ,ticks = LinearTicks(8)
    # ,color = :black
    ,label = "Pr=$(o.Pr), Ra=$(o.Ra)"
  )
end
map(outs[[8,10,13]]) do o
  cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
  lines!(xs, [evaluate!(cache, o.Nu, Gridap.Point([x,0])) for x in xs]
    ,ticks = LinearTicks(8)
    # ,color = :black
    ,label = "Pr=$(o.Pr), Ra=$(o.Ra)"
    ,linestyle=:dash
  )
end
axislegend(labelsize=10)
save("./nusselt_bot.pdf", f)
f
end
1+1
