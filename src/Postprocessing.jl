using Gridap
using Gridap.CellData
using Gridap.Arrays
using CairoMakie
using LineSearches
using Colors

using Random
Random.seed!(1234)

export plot_all_unitsquare, contourplot_unitsquare, contourplot_coannulus,
       cmap_cold_to_hot_paraview

function cmap_cold_to_hot_paraview()
  edges = [
    0.0, 0.45, 0.5, 0.55, 1.0
  ]
  colors = [
    RGBA(0, 1, 1, 1.0), RGBA(0, 0, 1, 1.0), RGBA(0, 0, 0.5019608, 1.0),
    RGBA(1, 0, 0, 1.0), RGBA(1, 1, 0, 1.0)
  ]
  return cgrad(colors, edges)
end

function contourplot_coannulus2(;
  fun::FEFunction=Gridap.interpolate(
    x->x[1]+x[2],
    TestFESpace(
      CartesianDiscreteModel((-2,2,-2,2),(10,10)),
      ReferenceFE(lagrangian, Float64, 2), conformity=:H1
    )
  ),
  name="out",
  ri = 2/3,
  ro = 5/3,
  eps = 0.01,
  n_plot=200,
  axisargs=(;),
  contourargs=(;),
  surfaceargs=(;),
  withsurface=true,
)

  rs = LinRange(ri+eps, ro-eps, n_plot)
  phis = LinRange(0, 2*pi, 2*n_plot-1)

  search_method = KDTreeSearch(num_nearest_vertices=5)
  fun_int = Interpolable(fun; searchmethod=search_method)
  cache = return_cache(fun_int, Gridap.Point(1.0, 1.0))
  function fun_help(x)
    return evaluate!(cache, fun_int, Gridap.Point(x))
  end
  xs = [r * cos(phi) for r in rs, phi in phis]
  ys = [r * sin(phi) for r in rs, phi in phis]
  zs = fun_help.(broadcast((x, y) -> (x, y), xs, ys))
  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    # title="stream function ψ",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(-ro, ro, -ro, ro);
    axisargs...
  )
  @show extrema(zs)
  if withsurface
    surface!(
      xs, ys, zs, shading = NoShading, colormap = :coolwarm, rasterize = true;
      surfaceargs...
    )
  end

  contour!(
    xs, ys, zs, levels = 10, labels=true, labelsize=15, color=:black;
    contourargs...
  )
  arc!(Point2f(0), ro, -π, π, color = :black, linewidth = 2)
  arc!(Point2f(0), ri, -π, π, color = :black, linewidth = 2)
  save("$(name).pdf", f)
  return f
end
contourplot_coannulus=contourplot_coannulus2


function contourplot_unitsquare2(;
  fun::FEFunction=Gridap.interpolate(
    x->x[1]+x[2],
    TestFESpace(
      CartesianDiscreteModel((0,1,0,1),(10,10)),
      ReferenceFE(lagrangian, Float64, 2), conformity=:H1
    )
  ),
  name="out",
  n_plot=100,
  axisargs=(;),
  contourargs=(;),
  surfaceargs=(;),
  withsurface=true,
)
  xs = LinRange(0, 1, n_plot)
  ys = LinRange(0, 1, n_plot)
  search_method = KDTreeSearch(num_nearest_vertices=10)

  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    # title="stream function ψ",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(0, 1, 0, 1);
    axisargs...
  )
  fun_int = Interpolable(fun; searchmethod=search_method)
  cache = return_cache(fun_int, Gridap.Point(0.0, 0.0))
  function fun_help(x)
    return evaluate!(cache, fun_int, Gridap.Point(x))
  end
  zs = [fun_help([x, y]) for x in xs, y in ys]
  @show extrema(zs)
  if withsurface
    surface!(
      xs, ys, zs, shading = NoShading, colormap = :coolwarm, rasterize = true;
      surfaceargs...
    )
  end
  contour!(
    xs, ys, zs, levels=10, labels=true, labelsize=15, color=:black,
    colormap = :coolwarm;
    contourargs...
  )
  save("$(name).pdf", f)
  return f
end
contourplot_unitsquare=contourplot_unitsquare2


"""
    plot_all_unitsquare4(psih, Pih, Th, uh, model, name, levels)

TBW
"""
function plot_all_unitsquare4(psih, Pih, Th, uh, model, name, levels)
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

  f = Figure(
    size=(500, 500), figure_padding=(0, 10, 0, 0)
  )
  Axis(
    f[1, 1],
    aspect=1,
    title="heat function Π",
    xlabel="x",
    ylabel="y",
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(6),
    yticks=LinearTicks(6),
    limits=(0, 1, 0, 1)
  )
  Pihi = Interpolable(Pih; searchmethod=search_method)
  cache_Pi = return_cache(Pihi, Gridap.Point(0.0, 0.0))
  function helper_Pi(x)
    return evaluate!(cache_Pi, Pihi, Gridap.Point(x))
    # TODO: Improve with using proper nodes and connectivity for plotting
  end
  zs = [helper_Pi([x, y]) for x in xs, y in ys]
  @info findmax(zs)

  contour!(
    xs, ys, zs, levels=levels.Pi, labels=true, labelsize=15, color=:black
  )
  save(joinpath(name, "heatfunction.pdf"), f)

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
  # writevtk(Ωₕ, joinpath(name, "nusselt.vtu"), cellfields=[
  #   "Nu" => Nu
  # ])

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

plot_all_unitsquare = plot_all_unitsquare4
