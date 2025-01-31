using GeMotion
using LineSearches: BackTracking, StrongWolfe, HagerZhang, MoreThuente
using Gridap
using Gridap.CellData
using Gridap.Arrays
using GridapGmsh
using CairoMakie

# TODO: Add tests

# Setup model
if haskey(ENV, "GITHUB_ACTIONS")
  model = GmshDiscreteModel(joinpath("co-annulus_unstructured_4.msh"))
  # model = GmshDiscreteModel(joinpath("co-annulus_structured_4.msh"))
else
  model = GmshDiscreteModel(joinpath("co-annulus_unstructured_4.msh"))
  # model = GmshDiscreteModel(joinpath("co-annulus_structured_5.msh"))
end


params = [
  [0.6, :newton, 6],
  [1.0, :trust_region, 8],
  [1.4, :newton, 8],
]

for p in params
  fname = "concentric-annulus" * "_$(p[1])"
  out = GeMotion.simulate(
    name = fname,
    Pr = 100,
    Ra = 1E4,
    n = p[1],
    model = model,
    nlsolver_opts = (;
      show_trace = true,
      method = p[2],
      linesearch = BackTracking(),
      ftol = 1E-8,
      xtol = 1E-50
    );
    (;
      T_diri_tags = ["inner", "outer"],
      T_diri_expressions = [1.0, 0.0]
    )...
  )

  # Postprocessing
  ri = 2/3
  ro = 5/3
  eps = 0.01
  rs = LinRange(ri+eps, ro-eps, 200)
  phis = LinRange(0, 2*pi, 2*200-1)

  search_method = KDTreeSearch(num_nearest_vertices=5)
  psihi = Interpolable(out.psih; searchmethod=search_method)
  cache_psi = return_cache(psihi, Gridap.Point(1.1, 0.1))
  function helper_psi(x)
    return evaluate!(cache_psi, psihi, Gridap.Point(x))
  end
  xs = [r * cos(phi) for r in rs, phi in phis]
  ys = [r * sin(phi) for r in rs, phi in phis]
  zs = helper_psi.(broadcast((x, y) -> (x, y), xs, ys))
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
    limits=(-ro, ro, -ro, ro)
  )
  surface!(xs, ys, zs, shading = NoShading, colormap = :coolwarm, rasterize = true)
  contour!(
    xs, ys, zs, levels = p[3] * 2, labels=true, labelsize=15, color=:black
  )
  arc!(Point2f(0), ro, -π, π, color = :black, linewidth = 2)
  arc!(Point2f(0), ri, -π, π, color = :black, linewidth = 2)
  display(f)
  save(joinpath(fname, "streamfunction.pdf"), f)


  search_method = KDTreeSearch(num_nearest_vertices=5)
  Thi = Interpolable(out.Th; searchmethod=search_method)
  cache_T = return_cache(Thi, Gridap.Point(0.0, 0.0))
  function helper_T(x)
    return evaluate!(cache_T, Thi, Gridap.Point(x))
  end
  zs = helper_T.(broadcast((x, y) -> (x, y), xs, ys))
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
    limits=(-ro, ro, -ro, ro)
  )
  # surface!(xs, ys, zs, shading = NoShading, colormap = :heat, rasterize = true)
  contour!(xs, ys, zs, levels=0.05:0.1:0.95, labels=true, labelsize=15, color=:black
  )
  contour!(xs, ys, zs, levels=0.1:0.1:0.90, colormap = :coolwarm
  )
  arc!(Point2f(0), ro, -π, π, color = :black, linewidth = 2)
  arc!(Point2f(0), ri, -π, π, color = :black, linewidth = 2)
  display(f)
  save(joinpath(fname, "temperature.pdf"), f)

end
