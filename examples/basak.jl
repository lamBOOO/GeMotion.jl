using GFluxx
using LineSearches: BackTracking, StrongWolfe
using Makie
using Gridap
using Gridap.Arrays  # return_cache

# solver settings
nlsolver_opts = (;
  show_trace=true,
  method=:newton,
  ftol=1E-8,
  xtol=1E-10
)

# BCs
uniform = (;
  T_diri_tags=[
    "botline", "leftline", "rightline", "botleftpoint", "botrightpoint"
  ],
  T_diri_expressions=[1.0, 0.0, 0.0, 0.5, 0.5]
)
wave = (;
  T_diri_tags=[
    "botline", "leftline", "rightline", "botleftpoint", "botrightpoint"
  ],
  T_diri_expressions=[x -> sin(pi * x[1]), 0.0, 0.0, 0.0, 0.0]
)

n_elems = 10  # 100 for paper
cases =
  [
    [0.7, 1E3, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([0.01, 0.05, 0.1, 0.15] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [0.7, 5 * 1E3, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([0.15, 0.5, 1, 1.3] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [0.7, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 5, 10, 13] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [0.1, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 4, 7, 9] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [1.0, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 5, 10, 14] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [10.0, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 5, 10, 14] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [0.015, 1E3, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([0.01, 0.05, 0.1, 0.15] |> x -> vcat(x, -x)), Sth=5, Sfl=5), uniform],
    [0.7, 1E3, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([0.01, 0.05, 0.1, 0.15] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave],
    [0.7, 5 * 1E3, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([0.15, 0.5, 1, 1.3] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave],
    [0.7, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 5, 10, 13] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave],
    [0.1, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 4, 7, 9] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave],
    [1.0, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 5, 10, 14] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave],
    [10.0, 1E5, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([1, 5, 10, 14] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave],
    [0.015, 1E3, 1.0, n_elems, (; nlsolver_opts..., linesearch=BackTracking()), (; T=[0.1 * i for i = 1:10], psi=([0.01, 0.05, 0.1, 0.15] |> x -> vcat(x, -x)), Sth=5, Sfl=5), wave]
  ]

# Run cases
outs = []
for case in cases
  out = GFluxx.simulate(Pr=case[1], Ra=case[2], n=case[3], n_elems=case[4], nlsolver_opts=case[5], levels=case[6]; case[7]...)
  push!(outs, out)
end


# Nusselt number plot
begin
  nn = 100
  xs = 0 |> x -> LinRange(0 + x, 1 - x, nn)
  update_theme!(
    palette=(color=Makie.wong_colors(), marker=[:circle, :xcross, :pentagon, :diamond],),
    Lines=(cycle=Cycle([:color, :marker], covary=true),)
  )
  f = Figure(
    # size = (500, 500)
    figure_padding=(0, 10, 0, 0)
  )
  ax = Axis(
    f[1, 1], aspect=300 / 180, title="local Nusselt number Nu, bottom wall", xlabel="x", ylabel="y", xminorticksvisible=true, yminorticksvisible=true, xticks=LinearTicks(5), yticks=LinearTicks(3), limits=((0.1, 0.9), (0, 15))
  )

  map(outs[[1, 3, 6]]) do o
    cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
    lines!(xs, [evaluate!(cache, o.Nu, Gridap.Point([x, 0])) for x in xs], ticks=LinearTicks(8)
      # ,color = :black
      , label="Pr=$(o.Pr), Ra=$(o.Ra)"
    )
  end
  map(outs[[8, 10, 13]]) do o
    cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
    lines!(xs, [evaluate!(cache, o.Nu, Gridap.Point([x, 0])) for x in xs], ticks=LinearTicks(8)
      # ,color = :black
      , label="Pr=$(o.Pr), Ra=$(o.Ra)", linestyle=:dash
    )
  end
  axislegend(labelsize=10)
  save("./nusselt_bot.pdf", f)
  f
end
1 + 1
