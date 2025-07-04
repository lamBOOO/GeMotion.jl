using GeMotion
using LineSearches: BackTracking, StrongWolfe, Static, HagerZhang, MoreThuente
using Makie
using Gridap
using Gridap.Arrays  # for "return_cache"
using CSV, DataFrames
using Printf

# Setup model
# (P3)-(L6)-(P4)
#  |         |
# (L7)      (L8)
#  |         |
# (P1)-(L5)-(P2)
model = CartesianDiscreteModel(
  (0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (40, 40) : (200, 200)
)

# solver settings
nlsolver_opts = (;
  show_trace=true,
  method=:trust_region,
  linesearch=BackTracking(),
  ftol=1E-8,
  xtol=1E-50
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

common_lvls = (; T=[0.1 * i for i = 1:10])

paramlist = [
  # Pr, Ra, n, psi_lvls, Pi_lvls, Sth_lvls, Sfl_lvls
  (0.7, 1E3, 1.0, [0.01; 0.05:0.05:0.15], 10, 5, 5, uniform)
  (0.7, 5 * 1E3, 1.0, [0.15, 0.5, 1, 1.3], 10, 5, 5, uniform)
  (0.7, 1E5, 1.0, [1, 5, 10, 13], 10, 5, 5, uniform)
  (0.1, 1E5, 1.0, [1, 4, 7, 9], 10, 5, 5, uniform)
  (1.0, 1E5, 1.0, [1, 5, 10, 14], 10, 5, 5, uniform)
  (10.0, 1E5, 1.0, [1, 5, 10, 14], 10, 5, 5, uniform)
  (
    0.015, 1E3, 1.0, [0.01; 0.05:0.05:0.15],
    vcat(1:1:4, 6, 8, 10:10:30, 0.01, 0.1, 0.5, 0.25),
    vcat(0.001, 0.003, 0.005, 0.01:0.01:0.05), 10, uniform
  )
  (0.7, 1E3, 1.0, [0.01, 0.05, 0.1, 0.15], 10, 5, 5, wave)
  (0.7, 5 * 1E3, 1.0, [0.15, 0.5, 1, 1.3], 10, 5, 5, wave)
  (0.7, 1E5, 1.0, [1, 5, 10, 13], 10, 5, 5, wave)
  (0.1, 1E5, 1.0, [1, 4, 7, 9], 10, 5, 5, wave)
  (1.0, 1E5, 1.0, [1, 5, 10, 14], 10, 5, 5, wave)
  (10.0, 1E5, 1.0, [1, 5, 10, 14], 10, 5, 5, wave)
  (0.015, 1E3, 1.0, [0.01, 0.05, 0.1, 0.15], 10, 5, 5, wave)
]
function mkcase(Pr, Ra, n, psi_vec, Pi_vec, Sth_vec, Sfl_vec, bcs)
  [
    Pr, Ra, n, model, nlsolver_opts,
    (;
      common_lvls...,
      psi=typeof(psi_vec) == Int ? psi_vec : vcat(psi_vec, -psi_vec),
      Pi=typeof(Pi_vec) == Int ? Pi_vec : vcat(Pi_vec, -Pi_vec),
      Sth=Sth_vec,
      Sfl=Sfl_vec
    ), bcs
  ]
end
cases = [
  mkcase(Pr, Ra, n, psi_vec, Pi_vec, Sth_vec, Sfl_vec, bcs)
  for (Pr, Ra, n, psi_vec, Pi_vec, Sth_vec, Sfl_vec, bcs) in paramlist
]

# Run cases
outs = []
outs2 = []
for (i, case) in enumerate(cases)
  out2 = GeMotion.simulate(
    name="$(i)_$(case[1])_$(case[2])", Pr=case[1], Ra=case[2], n=case[3],
    model=case[4], jac_scaling = 1, nlsolver_opts=case[5],
    nlsolver_custom_init_guess=[], nlsolver_init_guess_type=:zero,; case[7]...
  )
  out = GeMotion.plot_all_unitsquare(
    out2.psih, out2.Pih, out2.Th, out2.uh, model, "$(i)_$(case[1])_$(case[2])", case[6]
  )
  push!(outs2, out2)
  push!(outs, (; Pr=case[1], Ra=case[2], out...))
end

# Nusselt number plot bottom
refdatabot_fnames = [
  "basak_data/Uniform_Pr=0.7_Ra=1E3.txt",
  "basak_data/Uniform_Pr=0.7_Ra=1E5.txt",
  "basak_data/Uniform_Pr=10_Ra=1E5.txt",
  "basak_data/Non_Uniform_Pr=0.7_Ra=1E3.txt",
  "basak_data/Non_Uniform_Pr=0.7_Ra=1E5.txt",
  "basak_data/Non_Uniform_Pr=10_Ra=1E5.txt",
]
refdatabot = map(x -> CSV.read(x, DataFrame, header=false), refdatabot_fnames)
begin
  nn = 100
  xs = 0 |> x -> LinRange(0 + x, 1 - x, nn)
  update_theme!(
    palette=(
      color=Makie.wong_colors(), marker=[:circle, :xcross, :pentagon, :diamond],
    ),
    Lines=(cycle=Cycle([:color, :marker], covary=true),)
  )
  scale = 2.0
  f = Figure(
    size=scale .* (300, 180),
    # figure_padding= (0,0,0,0)
  )
  ax = Axis(
    f[1, 1],
    aspect=300 / 180,
    title="bottom wall",
    xlabel="position x",
    ylabel="Nusselt number Nu",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(5),
    yticks=LinearTicks(3),
    limits=((0.1, 0.9), (0, 15))
  )

  map(refdatabot) do data
    scatterlines!(
      data[!, 1],
      data[!, 2],
      marker=:circle,
      markersize=6,
      markercolor=(:black, 1.0),
      color=(:black, 1.0),
      # strokewidth=0.1,
    )
  end

  map(outs[[1, 3, 6]]) do o
    cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
    lines!(
      xs, [evaluate!(cache, o.Nu, Gridap.Point([x, 0])) for x in xs],
      label="Pr=$(o.Pr), Ra=$((@sprintf "%.0E" o.Ra))"
    )
  end
  map(outs[[8, 10, 13]]) do o
    cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
    lines!(
      xs, [evaluate!(cache, o.Nu, Gridap.Point([x, 0])) for x in xs],
      label="Pr=$(o.Pr), Ra=$((@sprintf "%.0E" o.Ra))", linestyle=:dash
    )
  end
  axislegend(
    labelsize=10,
    position=:ct,
    orientation=:horizontal,
    nbanks=3,
  )

  save("./nusselt_bot.pdf", f)
  f
end


# Nusselt number plot bottom
refdataside_fnames = [
  "basak_data/Side_Uniform_Pr=0.7_Ra=1E3.txt",
  "basak_data/Side_Uniform_Pr=0.7_Ra=1E5.txt",
  "basak_data/Side_Uniform_Pr=10_Ra=1E5.txt",
  "basak_data/Side_Non-Uniform_Pr=0.7_Ra=1E3.txt",
  "basak_data/Side_Non-Uniform_Pr=0.7_Ra=1E5.txt",
  "basak_data/Side_Non-Uniform_Pr=10_Ra=1E5.txt",
]
refdataside = map(x -> CSV.read(x, DataFrame, header=false), refdataside_fnames)
begin
  nn = 100
  xs = 0 |> x -> LinRange(0 + x, 1 - x, nn)
  update_theme!(
    palette=(
      color=Makie.wong_colors(), marker=[:circle, :xcross, :pentagon, :diamond],
    ),
    Lines=(cycle=Cycle([:color, :marker], covary=true),)
  )
  scale = 2.0
  f = Figure(
    size=scale .* (300, 180),
  )
  ax = Axis(
    f[1, 1],
    aspect=300 / 180,
    title="side wall",
    xlabel="position y",
    ylabel="Nusselt number Nu",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(5),
    yticks=LinearTicks(2),
    # hack to get same tick width:
    ytickformat=values -> ["  $(Int(value))" for value in values],
    limits=((0.1, 0.9), (0, 8))
  )

  map(refdataside) do data
    scatterlines!(
      data[!, 1],
      data[!, 2],
      marker=:circle,
      markersize=6,
      markercolor=(:black, 1.0),
      color=(:black, 1.0),
      # strokewidth=0.1,
    )
  end

  map(outs[[1, 3, 6]]) do o
    cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
    lines!(
      xs, -[evaluate!(cache, o.Nu, Gridap.Point([1, x])) for x in xs],
      label="Pr=$(o.Pr), Ra=$((@sprintf "%.0E" o.Ra))"
    )
  end
  map(outs[[8, 10, 13]]) do o
    cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
    lines!(
      xs, -[evaluate!(cache, o.Nu, Gridap.Point([1, x])) for x in xs],
      label="Pr=$(o.Pr), Ra=$((@sprintf "%.0E" o.Ra))", linestyle=:dash
    )
  end
  axislegend(
    labelsize=10,
    position=:ct,
    orientation=:horizontal,
    nbanks=3,
  )

  save("./nusselt_side.pdf", f)
  f
end

@info "relative errors bot"
for i = 1:6
  # begin; i = 3
  i_sim = [1, 3, 6, 8, 10, 13][i]
  otmp = outs[i_sim]
  cache_tmp = return_cache(otmp.Nu, Gridap.Point(0.0, 0.0))
  results_sim = [
    evaluate!(cache_tmp, otmp.Nu, Gridap.Point([x, 0]))
    for x in refdatabot[i][!, 1]
  ]
  diff = abs.(results_sim - refdatabot[i][!, 2])
  l2norm_ref = sqrt(sum(
    (refdatabot[i][!, 1] |> x -> [x[i+1] - x[i] for i in 1:length(x)-1])
    .* refdatabot[i][!, 2][1:end-1] .^ 2
  ))
  l2norm_diff = sqrt(sum(
    (refdatabot[i][!, 1] |> x -> [x[i+1] - x[i] for i in 1:length(x)-1])
    .* diff[1:end-1] .^ 2
  ))
  rel_l2err = l2norm_diff / l2norm_ref
  @show rel_l2err
  @assert rel_l2err < 0.0312
  lines(refdatabot[i][!, 1], diff, color=:black)
end

@info "relative errors side"
for i = 1:6
  # begin; i = 1
  i_sim = [1, 3, 6, 8, 10, 13][i]
  otmp = outs[i_sim]
  cache_tmp = return_cache(otmp.Nu, Gridap.Point(0.0, 0.0))
  results_sim = -[
    evaluate!(cache_tmp, otmp.Nu, Gridap.Point([1, x]))
    for x in refdataside[i][!, 1]
  ]
  diff = abs.(results_sim - refdataside[i][!, 2])
  l2norm_ref = sqrt(sum(
    (refdataside[i][!, 1] |> x -> [x[i+1] - x[i] for i in 1:length(x)-1])
    .* refdataside[i][!, 2][1:end-1] .^ 2
  ))
  l2norm_diff = sqrt(sum(
    (refdataside[i][!, 1] |> x -> [x[i+1] - x[i] for i in 1:length(x)-1])
    .* diff[1:end-1] .^ 2
  ))
  rel_l2err = l2norm_diff / l2norm_ref
  @show rel_l2err
  @assert rel_l2err < 0.06
  lines(refdataside[i][!, 1], diff, color=:black)
end
