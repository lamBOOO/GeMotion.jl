using GeMotion
using LineSearches: BackTracking, StrongWolfe, HagerZhang, MoreThuente, Static, InitialStatic
using Gridap
using Gridap.CellData
using Gridap.Arrays
using GridapGmsh
using CairoMakie
using Colors
using FileIO

# 1)
model_square = CartesianDiscreteModel(
(0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (40, 40) : (20, 20)
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
  Pr=100,
  Ra=1e4,
  n=0.6,
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
    # alphaguess=InitialStatic(alpha=1.0),
    linesearch=BackTracking(
      # order=3,
      # c_1=0.0,
      # ρ_hi=0.99,
      # ρ_lo=0.01,
      # iterations=2
    ),
    # linesearch=HagerZhang(
    #   display=1
    # ),
    # linesearch=Static(scaled=true),
    # linesearch=HagerZhang(),
    ftol=1E-8,
    xtol=1E-50,
    iterations=50,
    # extended_trace=true
  );
  nlsolver_custom_init_guess=[],
  # nlsolver_custom_init_guess=out1.result.free_values,
  # nlsolver_custom_init_guess=vcat(
  #   zeros(size(out1.result[1].free_values)),
  #   out1.result[2].free_values,
  #   out1.result[3].free_values,
  # ),
  nlsolver_init_guess_type=:zero,
  # nlsolver_init_guess_type=:custom,
  wave...
)

out_post = GeMotion.plot_all_unitsquare(
  out1.psih, out1.Th, out1.uh, model_square, name, (T=[0.1 * i for i = 1:10], psi=[0.01,0.05,0.1,0.5,1.5,3.5,5.5,6]|>x->vcat(-x,x), Sth=[0.1,0.2,0.3,0.4,0.5,1.0,2,5], Sfl=vcat([0.1,1,10,50,100]))
)
1+1

