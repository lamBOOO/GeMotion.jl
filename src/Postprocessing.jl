using Gridap
using Gridap.CellData
using Gridap.Arrays
using CairoMakie
using LineSearches

using Random
Random.seed!(1234)

export plot_all_unitsquare

"""
    plot_all_unitsquare2(psih, Th, uh, model, name, levels)

TBW
"""
function plot_all_unitsquare2(psih, Th, uh, model, name, levels)
  # TODO: Make named args

  Ωₕ = Triangulation(model)
  btrian = BoundaryTriangulation(model)

  # plotting
  n_plot = 100
  xs = LinRange(0, 1, n_plot)
  ys = LinRange(0, 1, n_plot)
  search_method = KDTreeSearch(num_nearest_vertices=10)

  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    title="stream function ψ",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(0, 1, 0, 1)
  )
  psihi = Interpolable(psih; searchmethod=search_method)
  cache_psi = return_cache(psihi, Gridap.Point(0.0, 0.0))
  function helper_psi(x)
    return evaluate!(cache_psi, psihi, Gridap.Point(x))
    # TODO: Improve with using proper nodes and connectivity for plotting
  end
  zs = [helper_psi([x, y]) for x in xs, y in ys]
  @info findmax(zs)

  contour!(
    xs, ys, zs, levels=levels.psi, labels=true, labelsize=15, color=:black
  )
  save(joinpath(name, "streamfunction.pdf"), f)

  Thi = Interpolable(Th; searchmethod=search_method)
  cache_T = return_cache(Thi, Gridap.Point(0.0, 0.0))
  function helper_T(x)
    return evaluate!(cache_T, Thi, Gridap.Point(x))
  end
  zs = [helper_T([x, y]) for x in xs, y in ys]
  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    title="temperature T",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(0, 1, 0, 1)
  )
  contour!(xs, ys, zs, levels=levels.T, labels=true, labelsize=15, color=:black
  )
  save(joinpath(name, "temperature.pdf"), f)

  # Nusselt number
  # TODO: Add left wall
  nb = get_normal_vector(btrian)
  Nu = Interpolable(
    (∇(Th) ⋅ nb); searchmethod=KDTreeSearch(num_nearest_vertices=5)
  )

  # entropies
  Sth = interpolate(
    ∇(Th) |> x -> inner(x, x),
    TestFESpace(model, ReferenceFE(lagrangian, Float64, 2), conformity=:H1)
  )
  Sfl = interpolate(∇(uh) |> x -> 1E-4 * (
      2 * (inner(x .* x, TensorValue(1, 0, 0, 1)))
      +
      ((inner(x, TensorValue(0, 1, 1, 0))) |> x -> x * x)
    ), TestFESpace(model, ReferenceFE(lagrangian, Float64, 2), conformity=:H1))
  writevtk(Ωₕ, joinpath(name, "entropy.vtu"), cellfields=[
    "Sth" => Sth, "Sfl" => Sfl
  ])

  xs = LinRange(0, 1, n_plot)
  ys = LinRange(0, 1, n_plot)

  cache_Sth = return_cache(Sth, Gridap.Point(0.0, 0.0))
  cache_Sfl = return_cache(Sfl, Gridap.Point(0.0, 0.0))
  z_Sfl = [evaluate!(cache_Sfl, Sfl, Gridap.Point([x, y])) for x in xs, y in ys]
  z_Sth = [evaluate!(cache_Sth, Sth, Gridap.Point([x, y])) for x in xs, y in ys]

  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    title="local heat entropy S_θ",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(0, 1, 0, 1)
  )
  contour!(xs, ys, z_Sth, levels=levels.Sth
    # ,levels=vcat(1:1:7, 10:5:30, 0.01,0.1,0.5,0.25)
    # ,levels=10
    , labels=true, labelsize=15, color=:black
  )
  save(joinpath(name, "Sth.pdf"), f)

  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    title="local fluid entropy S_ψ",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(0, 1, 0, 1)
  )
  function formatter(level::Real)::String
    lev_short = round(level; digits=3)
    string(isinteger(lev_short) ? round(Int, lev_short) : lev_short)
  end
  contour!(xs, ys, z_Sfl, levels=levels.Sfl
    # ,levels=vcat(0.005:0.005:0.05)
    # ,levels=10
    , labels=true, labelsize=15, labelformatter=formatter, color=:black
  )
  save(joinpath(name, "Sfl.pdf"), f)

  return (
    Nu=Nu,
    Sth=Sth,
    Sfl=Sfl,
    btrian=btrian,
    model=model,
    Ωₕ=Ωₕ,
  )
end

plot_all_unitsquare = plot_all_unitsquare2
