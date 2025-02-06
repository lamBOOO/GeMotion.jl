using GeMotion
using LineSearches: BackTracking, StrongWolfe, HagerZhang, MoreThuente
using Gridap
using Gridap.CellData
using Gridap.Arrays
using GridapGmsh
using CairoMakie
using Colors
using FileIO

# Setup model
if haskey(ENV, "GITHUB_ACTIONS")
  model = GmshDiscreteModel(joinpath("meshes/co-annulus_unstructured_4.msh"))
  # model = GmshDiscreteModel(joinpath("co-annulus_structured_3.msh"))
else
  model = GmshDiscreteModel(joinpath("meshes/co-annulus_unstructured_4.msh"))
  # model = GmshDiscreteModel(joinpath("co-annulus_structured_4.msh"))
end

name = "out"

out = GeMotion.simulate(
  name=name,
  Pr=0.706,
  Ra=4.7e4,
  n=1.0,
  model=model,
  nlsolver_opts=(;
    show_trace=true,
    method=:newton,
    linesearch=BackTracking(),
    ftol=1E-8,
    xtol=1E-50
  );
  (;
    T_diri_tags=["inner", "outer"],
    T_diri_expressions=[1.0, 0.0],
    V_diri_tags=["all"]
  )...
)
1 + 1

n_plot = 200
ri = 5 / 8
ro = 13 / 8
eps = 0.001
rs = LinRange(ri + eps, ro - eps, n_plot)
phis = LinRange(0, 2 * pi, 2 * n_plot - 1)
xs = [r * cos(phi) for r in rs, phi in phis]
ys = [r * sin(phi) for r in rs, phi in phis]

search_method = KDTreeSearch(num_nearest_vertices=5)
Thi = Interpolable(out.Th; searchmethod=search_method)
cache_T = return_cache(Thi, Gridap.Point(0.0, 0.0))
function helper_T(x)
  return evaluate!(cache_T, Thi, Gridap.Point(x))
end
zs = helper_T.(broadcast((x, y) -> (x, y), xs, ys))

begin
  f = Figure(
    size=(2 * 500 - 65, 500), figure_padding=(-10, 0, 0, 0), rowgap=0, colgap=0
  )
  ax = Axis(
    f[1, 1],
    aspect=1,
    title="Simulation",
    titlesize=20,
    xlabel="x",
    ylabel="y",
    xlabelsize=20,
    ylabelsize=20,
    xminorticksvisible=true,
    yminorticksvisible=true,
    xticks=LinearTicks(8),
    yticks=LinearTicks(8),
    limits=(-ro, ro, -ro, ro),
  )

  edges = [
    0.0, 0.035, 0.08, 0.2, 0.25, 0.34, 0.38, 0.44, 0.48, 0.54, 0.58, 0.65,
    0.74, 0.83, 0.88, 0.95, 1.0
  ]
  colors = [
    isodd(i) ? RGBA(1, 1, 1, 0.0) : RGBA(0, 0, 0, 0.05)
    for i in 1:(length(edges)-1)
  ]
  cmap = cgrad(colors, edges; categorical=true)

  surface!(
    xs, ys, zs,
    shading=NoShading,
    rasterize=true,
    colormap=:coolwarm,
    colorrange=(0, 1),
  )
  surface!(
    xs, ys, zs,
    shading=NoShading,
    rasterize=true,
    colormap=cmap,
    colorrange=(0, 1),
    transparency=true,
    # alpha = 0.2,
  )
  contour!(
    xs, ys, zs,
    levels=edges,
    labels=true, labelsize=15, color=:black
  )
  # contour!(xs, ys, zs, levels=0.1:0.1:0.90, colormap = :coolwarm
  # )
  arc!(Point2f(0), ro, -π, π, color=:red, linewidth=2)
  arc!(Point2f(0), ri, -π, π, color=:red, linewidth=2)
  arc!(Point2f(0), 0.02, -π, π, color=:red, linewidth=5)

  img = load("experiment/ex1.png")
  ax2 = Axis(f[1, 2];
    xticksvisible=false,
    yticksvisible=false,
    xlabelvisible=false,
    ylabelvisible=false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    aspect=DataAspect(),
    title="Experiment [Kuehn/Goldstein 1976]",
    titlesize=20
  )
  img = img[15:end-15, 60:end-70]
  image!(ax2, -ro .. ro, -ro .. ro, rotr90(img))
  @info size(img)
  @assert size(img)[1] == size(img)[2]
  arc!(ax2, Point2f(0), ro, -π, π, color=:red, linewidth=2)
  arc!(ax2, Point2f(0), ri, -π, π, color=:red, linewidth=2)
  arc!(ax2, Point2f(0), 0.02, -π, π, color=:red, linewidth=5)

  display(f)
  save(joinpath(name, "temperature_comparison.pdf"), f)
  save(joinpath(name, "temperature_comparison.png"), f)

  @assert abs(
    [helper_T([r, 0]) for r in rs] |> x -> sum(x) / length(x) - 0.33
  ) < 0.01 "mean of T(x,0) is close to 0.33"
end


