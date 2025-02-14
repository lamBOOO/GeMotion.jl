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
(0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (40, 40) : (120, 120)
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


paramlist = [
  # n, model, bcs
  (0.6, model_square, uniform)
  (1.0, model_square, uniform)
  (1.4, model_square, uniform)
  (0.6, model_square, wave)
  (1.0, model_square, wave)
  (1.4, model_square, wave)
]
function mkcase(n, model, bcs)
  (;n, model, bcs)
end
cases = [
  mkcase(n, model, bcs)
  for (n, model, bcs) in paramlist
]

# outs1 = []
out = nothing
for (i, case) in enumerate(cases)
  name = "square_$(i)_$(case.n)"
  out = GeMotion.simulate(
    name=name,
    Pr=100,
    Ra=1e4,
    n=case.n,
    model=case.model,
    jac_scaling = 1,
    nlsolver_opts=(;
      show_trace=true,
      method=:newton,
      linesearch=BackTracking(),
      ftol=1E-8,
      xtol=1E-50,
      iterations=50,
    );
    nlsolver_custom_init_guess=[],
    nlsolver_init_guess_type=:zero,
    case.bcs...
  )
  push!(outs1, out)
  GeMotion.contourplot_unitsquare(;
    fun=out.psih,
    name=name*"/psi",
    contourargs=(;levels=[2.0^i for i in -1:1:3]|>x->vcat(-x,x)),
    surfaceargs=(;
      colormap=:coolwarm
    )
  ) |> display
  GeMotion.contourplot_unitsquare(;
    fun=out.Th,
    name=name*"/T",
    contourargs=(;levels=0:0.1:1),
    surfaceargs=(;
      colormap=cmap_cold_to_hot_paraview()
    )
  ) |> display
  GeMotion.contourplot_unitsquare(;
    fun=out.Sth,
    name=name*"/Sth",
    contourargs=(;levels=[4.0^i for i in -1:1:3]),
    surfaceargs=(;
    colormap=:Reds,
      # colorscale=x->min(x,100)
    )
  ) |> display
  GeMotion.contourplot_unitsquare(;
    fun=out.Sfl,
    name=name*"/Sfl",
    contourargs=(;levels=[4.0^i for i in -2:1:3]),
    surfaceargs=(;
    colormap=:Blues,
      # colorrange=(0,100)
    )
  ) |> display
end

begin
  scale = 1.5
  f = Figure(
    size=scale .* (600, 250),
    figure_padding= (0,10,0,0)
  )
  nplot_Nu = 100
  xs_Nu = LinRange(0.0, 1.0, 2*nplot_Nu)[1:end]

  map(enumerate([1:3,4:6])) do (i_indices, indices)

    ax = Axis(
      f[1, i_indices],
      title=i_indices==1 ? "inner wall, uniform heating" : "inner wall, non-uniform heating",
      xlabel="angle ϕ",
      ylabel="Nusselt number Nu",
      xminorticksvisible=true,
      yminorticksvisible=true,
      xticks=LinearTicks(5),
      yticks=LinearTicks(5),
      limits=((0.0, 1.0), (0,20))
    )

    map(outs1[indices]) do o
      @info "hi"
      nb = get_normal_vector(o.btrian)
      NUU = Interpolable(
        # get (-dT/dr) with a transformation to polar coordinates
        (- ∇(o.Th) ⋅ (x->VectorValue([0.0,1.0])));
          searchmethod=KDTreeSearch(num_nearest_vertices=500)
      )
      cache = return_cache(NUU, Gridap.Point(0.0, 0.0))
      Nu_inner = [
        evaluate!(cache, NUU, Gridap.Point([x, 0.0]))
        for x in xs_Nu
      ]
      lines!(
        xs_Nu,
        Nu_inner,
        label="n = $(o.n)",
      )
    end
    axislegend(
      labelsize=10,
      position=:lt,
      orientation=:horizontal,
      nbanks=1,
    )
  end
  f |> display
  save("nusselt_number_square.pdf", f)
end



# 2)
if haskey(ENV, "GITHUB_ACTIONS")
  model_annulus = GmshDiscreteModel(
    joinpath("../meshes/2.5/co-annulus_unstructured_4.msh")
  )
  # model_annulus = GmshDiscreteModel(
  #   joinpath("../meshes/2.5/co-annulus_structured_2.msh")
  # )
else
  model_annulus = GmshDiscreteModel(
    joinpath("../meshes/2.5/co-annulus_unstructured_6.msh")
  )
  # model_annulus = GmshDiscreteModel(
  #   joinpath("../meshes/2.5/co-annulus_structured_6.msh")
  # )
end
uniform_annulus = (;
  T_diri_tags=["inner", "outer"],
  T_diri_expressions=[1.0, 0.0],
  V_diri_tags=["all"]
)
wave_annulus = (;
  T_diri_tags=["inner", "outer"],
  T_diri_expressions=[x->0.5*(sin(atan(x[2],x[1]))+1), 0.0],
  V_diri_tags=["all"]
)

paramlist = [
  # n, model, bcs, Sfl_lvls
  (0.6, model_annulus, uniform_annulus, [10.0^i for i=-2:1:2])
  (1.0, model_annulus, uniform_annulus, [4.0^i for i in -0:1:3])
  (1.4, model_annulus, uniform_annulus, [4.0^i for i in -0:1:3])
  (0.6, model_annulus, wave_annulus, [10.0^i for i=-2:1:2])
  (1.0, model_annulus, wave_annulus, [4.0^i for i in -2:1:3])
  (1.4, model_annulus, wave_annulus, [4.0^i for i in -2:1:3])
]
function mkcase(n, model, bcs, Sfl_lvls)
  (;n, model, bcs, Sfl_lvls)
end
cases = [
  mkcase(n, model, bcs, Sfl_lvls)
  for (n, model, bcs, Sfl_lvls) in paramlist
]

# outs2 = []
for (i, case) in enumerate(cases)
  name = "annulus_$(i)_$(case.n)"

  # initialize with a slightly increased n
  # this avoids "bad" solution for uniform n=1.0
  out = GeMotion.simulate(
    name=name,
    Pr=100,
    Ra=1e4,
    n=case.n+0.1,
    model=case.model,
    jac_scaling = 1,
    nlsolver_opts=(;
      show_trace=true,
      method=:newton,
      linesearch=BackTracking(),
      ftol=1E-8,
      xtol=1E-50,
      iterations=200,
    );
    nlsolver_custom_init_guess=[],
    nlsolver_init_guess_type=:zero,
    case.bcs...
  )

  out = GeMotion.simulate(
    name=name,
    Pr=100,
    Ra=1e4,
    n=case.n,
    model=case.model,
    jac_scaling = 1,
    nlsolver_opts=(;
      show_trace=true,
      method=:newton,
      linesearch=BackTracking(),
      ftol=1E-8,
      xtol=1E-50,
      iterations=200,
    );
    nlsolver_custom_init_guess=out.result.free_values,
    nlsolver_init_guess_type=:custom,
    case.bcs...
  )
  push!(outs2, out)
  opts=(;ri=2/3,ro=5/3,eps=0.01,n_plot=100)
  GeMotion.contourplot_coannulus(
    ;opts...,
    fun=out.psih,
    name=name*"/psi",
    contourargs=(;levels=[2.0^i for i in -1:1:3]|>x->vcat(-x,x)),
    surfaceargs=(;colormap=:coolwarm)
  ) |> display
  GeMotion.contourplot_coannulus(
    ;opts...,
    fun=out.Th,
    name=name*"/T",
    contourargs=(;levels=0:0.1:1),
    surfaceargs=(;colormap=cmap_cold_to_hot_paraview())
  ) |> display
  GeMotion.contourplot_coannulus(
    ;opts...,
    fun=out.Sth,
    name=name*"/Sth",
    contourargs=(;levels=[4.0^i for i in -1:1:3]),
    surfaceargs=(;
      colormap=:Reds,
    )
  ) |> display
  GeMotion.contourplot_coannulus(
    ;opts...,
    fun=out.Sfl,
    name=name*"/Sfl",
    contourargs=(;levels=case.Sfl_lvls),
    surfaceargs=(;colormap=:Blues)
  ) |> display
end

begin
  scale = 1.5
  f = Figure(
    size=scale .* (600, 250),
    figure_padding= (0,10,0,0)
  )
  nplot_Nu = 100
  phis_Nu = LinRange(0, 2*pi, 2*nplot_Nu)[1:end]

  map(enumerate([1:3,4:6])) do (i_indices, indices)

    ax = Axis(
      f[1, i_indices],
      title=i_indices==1 ? "inner wall, uniform heating" : "inner wall, non-uniform heating",
      xlabel="angle ϕ",
      ylabel="Nusselt number Nu",
      xminorticksvisible=true,
      yminorticksvisible=true,
      xticks = (0:1/2*pi:2*pi, ["0", "π/2", "π", "3π/2", "2π"]),
      yticks=LinearTicks(5),
      limits=((0.0, 2*pi), (-0.5,10.5))
    )

    map(outs2[indices]) do o
      nb = get_normal_vector(o.btrian)
      Nu = Interpolable(
        # get (-dT/dr) with a transformation to polar coordinates
        (- ∇(o.Th) ⋅ (x->VectorValue([x[1],x[2]] / sqrt(x[1]^2 + x[2]^2))));
          searchmethod=KDTreeSearch(num_nearest_vertices=50)
      )
      cache = return_cache(Nu, Gridap.Point(0.0, 0.0))
      dist = 0.0
      ri = 2/3
      Nu_inner = [
        evaluate!(cache, Nu, Gridap.Point([(ri+dist)*cos(p), (ri+dist)*sin(p)]))
        for p in phis_Nu
      ]
      lines!(
        phis_Nu,
        Nu_inner,
        label="n = $(o.n)",
      )
    end
    axislegend(
      labelsize=10,
      position=:lt,
      orientation=:horizontal,
      nbanks=1,
    )
  end
  f |> display
  save("nusselt_number_annulus.pdf", f)
end
