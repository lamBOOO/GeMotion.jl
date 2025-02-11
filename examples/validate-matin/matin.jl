using GeMotion
using LineSearches: BackTracking, StrongWolfe, HagerZhang, MoreThuente
using Gridap
using Gridap.CellData
using Gridap.Arrays
using GridapGmsh
using CairoMakie
using Colors

# TODO: Add tests

# Setup model
if haskey(ENV, "GITHUB_ACTIONS")
  model = GmshDiscreteModel(
    joinpath("../meshes/2.5/co-annulus_unstructured_4.msh")
  )
  # model = GmshDiscreteModel(
  #   joinpath("../meshes/2.5/co-annulus_structured_4.msh")
  # )
else
  model = GmshDiscreteModel(
    joinpath("../meshes/2.5/co-annulus_unstructured_4.msh")
  )
  # model = GmshDiscreteModel(
  #   joinpath("../meshes/2.5/co-annulus_structured_4.msh")
  # )
end


params = [
  [0.6, :newton, 6],
  [0.7, :newton, 8],
  [0.8, :newton, 8],
  [0.9, :newton, 8],
  [1.0, :newton, 8],
  [1.1, :newton, 8],
  [1.2, :newton, 8],
  [1.3, :newton, 8],
  [1.4, :newton, 8],
]

nusselts = []
for (i, p) in enumerate(params)
  fname = "$(i)_$(p[1])_$(string(p[2]))"
  out = GeMotion.simulate(
    name = fname,
    Pr = 100,
    Ra = 1E4,
    n = p[1],
    model = model,
    jac_scaling = 1e-6,
    nlsolver_opts = (;
      show_trace = true,
      method = p[2],
      linesearch = BackTracking(),
      ftol = 1E-10,
      xtol = 1E-20
    );
    nlsolver_custom_init_guess=[],
    nlsolver_init_guess_type=:zero,
    (;
      T_diri_tags = ["inner", "outer"],
      T_diri_expressions = [1.0, 0.0],
      V_diri_tags = ["all"]
    )...
  )

  # Postprocessing
  edges = [
    0.0, 0.45, 0.5, 0.55, 1.0
  ]
  colors = [
    RGBA(0, 1, 1, 1.0), RGBA(0, 0, 1, 1.0), RGBA(0, 0, 0.5019608, 1.0),
    RGBA(1, 0, 0, 1.0), RGBA(1, 1, 0, 1.0)
  ]

  cmap_cold_to_hot_paraview = cgrad(colors, edges)
  ri = 2/3
  ro = 5/3
  eps = haskey(ENV, "GITHUB_ACTIONS") ? 0.1 : 0.001
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
  surface!(
    xs, ys, zs, shading = NoShading, colormap = :coolwarm, rasterize = true
  )
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
  surface!(
    xs, ys, zs, shading = NoShading, colormap = cmap_cold_to_hot_paraview,
    rasterize = true
  )
  contour!(
    xs, ys, zs, levels=0.05:0.1:0.95, labels=true, labelsize=15, color=:black
  )
  # contour!(xs, ys, zs, levels=0.1:0.1:0.90, colormap = :coolwarm)
  arc!(Point2f(0), ro, -π, π, color = :black, linewidth = 2)
  arc!(Point2f(0), ri, -π, π, color = :black, linewidth = 2)
  display(f)
  save(joinpath(fname, "temperature.pdf"), f)


  # mean Nusselt number over inner cylinder
  nb = get_normal_vector(out.btrian)
  Nu = Interpolable(
    # get (-dT/dr) with a transformation to polar coordinates
    (- ∇(out.Th) ⋅ (x->VectorValue([x[1],x[2]] / sqrt(x[1]^2 + x[2]^2))));
     searchmethod=KDTreeSearch(num_nearest_vertices=5)
  )
  cache = return_cache(Nu, Gridap.Point(0.0, 0.0))
  phis_Nu = LinRange(0, 2*pi, 2*1000-1)[1:end-1]
  nussel_number = [
    evaluate!(cache, Nu, Gridap.Point([(2/3+0.0)*cos(p), (2/3+0.0)*sin(p)]))
    for p in phis_Nu
  ] |> x -> sum(x)/length(x)
  println("Nusselt number at inner cylinder:", nussel_number)
      push!(nusselts, nussel_number)

end


# begin
#   f = Figure(
#     size=(250, 215), figure_padding=(0, 10, 0, 0),
#   )
#   Axis(
#     f[1, 1],
#     aspect=250/215,
#     xticks=0.6:0.2:1.4,
#     yticks=1.0:1.5:8.5,
#     limits = (0.6, 1.4, 1.0, 8.5)
#   )
#   scatterlines!(
#     collect(0.6:0.1:1.4), nusselts
#   )
#   display(f)
# end
