using GeMotion
using LineSearches: BackTracking, StrongWolfe, HagerZhang, MoreThuente
using Gridap
using Gridap.CellData
using Gridap.Arrays
using GridapGmsh
using CairoMakie
using Colors
using FileIO

# 1)
model_square = CartesianDiscreteModel(
  (0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (40, 40) : (80, 80)
)
# BCs
uniform = (;
  T_diri_tags=[5, 7, 8, 1, 2],
  T_diri_expressions=[1.0, 0.0, 0.0, 0.5, 0.5],
  V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
)
wave = (;
  T_diri_tags=[5, 7, 8, 1, 2],
  T_diri_expressions=[x -> sin(pi * x[1]), 0.0, 0.0, 0.0, 0.0],
  V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8]
)

# 2)
# if haskey(ENV, "GITHUB_ACTIONS")
#   model_annulus = GmshDiscreteModel(
#     joinpath("../meshes/2.6/co-annulus_unstructured_4.msh")
#   )
#   # model_annulus = GmshDiscreteModel(
#   #   joinpath("../meshes/2.6/co-annulus_structured_3.msh")
#   # )
# else
#   # model_annulus = GmshDiscreteModel(
#   #   joinpath("../meshes/2.6/co-annulus_unstructured_2.msh")
#   # )
#   model_annulus = GmshDiscreteModel(
#     joinpath("../meshes/2.6/co-annulus_structured_2.msh")
#   )
# end

name = "con-annulus"

out1 = GeMotion.simulate(
  name=name,
  Pr=0.026,
  Ra=1e5,
  n=1.0,
  model=model_square,
  jac_scaling = 1,
  nlsolver_opts=(;
    show_trace=true,
    # method=:trust_region,
    # method=:anderson,
    # beta=0.01,
    # factor=0.1,
    # autoscale=false,
    method=:newton,
    linesearch=BackTracking(),
    # linesearch=HagerZhang(),
    ftol=1E-8,
    xtol=1E-50,
    iterations=50,
    # extended_trace=true
  );
  nlsolver_custom_init_guess=[],
  nlsolver_init_guess_type=:zero,1
  wave...
)

out_post = GeMotion.plot_all_unitsquare(
  out1.psih, out1.Th, out1.uh, model_square, name, (T=[0.1 * i for i = 1:10], psi=[0.01,0.05,0.1,0.5,1.5,3.5,5.5,6]|>x->vcat(-x,x), Sth=10, Sfl=10)
)
1+1






# n_plot = 200
# ri = 5 / 8
# ro = 13 / 8
# eps = 0.01
# rs = LinRange(ri + eps, ro - eps, n_plot)
# phis = LinRange(0, 2 * pi, 2 * n_plot - 1)
# xs = [r * cos(phi) for r in rs, phi in phis]
# ys = [r * sin(phi) for r in rs, phi in phis]

# search_method = KDTreeSearch(num_nearest_vertices=5)
# Thi = Interpolable(out.Th; searchmethod=search_method)
# cache_T = return_cache(Thi, Gridap.Point(0.0, 0.0))
# function helper_T(x)
#   return evaluate!(cache_T, Thi, Gridap.Point(x))
# end
# zs = helper_T.(broadcast((x, y) -> (x, y), xs, ys))

# begin
#   search_method = KDTreeSearch(num_nearest_vertices=5)
#   psihi = Interpolable(out.psih; searchmethod=search_method)
#   cache_psi = return_cache(psihi, Gridap.Point(1.1, 0.1))
#   function helper_psi(x)
#     return evaluate!(cache_psi, psihi, Gridap.Point(x))
#   end
#   xs = [r * cos(phi) for r in rs, phi in phis]
#   ys = [r * sin(phi) for r in rs, phi in phis]
#   zs = helper_psi.(broadcast((x, y) -> (x, y), xs, ys))
#   f = Figure(
#     size=(500, 500), figure_padding=(0, 10, 0, 0)
#   )
#   Axis(
#     f[1, 1],
#     aspect=1,
#     title="stream function ψ",
#     xlabel="x",
#     ylabel="y",
#     xminorticksvisible=true,
#     yminorticksvisible=true,
#     xticks=LinearTicks(6),
#     yticks=LinearTicks(6),
#     limits=(-ro, ro, -ro, ro)
#   )
#   surface!(
#     xs, ys, zs, shading = NoShading, colormap = :coolwarm, rasterize = true
#   )
#   contour!(
#     xs, ys, zs, levels = 10, labels=true, labelsize=15, color=:black
#   )
#   arc!(Point2f(0), ro, -π, π, color = :black, linewidth = 2)
#   arc!(Point2f(0), ri, -π, π, color = :black, linewidth = 2)
#   display(f)
# end


