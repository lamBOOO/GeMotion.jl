using GeMotion
using LineSearches: BackTracking, StrongWolfe, HagerZhang, MoreThuente, Static, InitialStatic
using Gridap
using Gridap.CellData
using Gridap.Arrays
using GridapGmsh
using CairoMakie
using Colors
using FileIO
using CSV, DataFrames

do_conv_study = haskey(ENV, "GITHUB_ACTIONS") ? false : true
do_study = haskey(ENV, "GITHUB_ACTIONS") ? true : false


# 1)
# Setup model
# (P3)-(L6)-(P4)
#  |         |
# (L7)      (L8)
#  |         |
# (P1)-(L5)-(P2)
model_square = CartesianDiscreteModel(
  (0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (40, 40) :
  round.(Int, floor.((61, 61).*sqrt(2)^(i-1)))
)
model_squares = [CartesianDiscreteModel(
  (0, 1, 0, 1), haskey(ENV, "GITHUB_ACTIONS") ? (40, 40) :
  round.(Int, floor.((61, 61).*sqrt(2)^(i-1))) |> x-> x .+ isodd.(x)  # even num
) for i in 1:5]
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

paramlist = [
  # n, model, bcs
  do_study ? [
    (0.6, model_square, uniform)
    (1.0, model_square, uniform)
    (1.4, model_square, uniform)
    (0.6, model_square, wave)
    (1.0, model_square, wave)
    (1.4, model_square, wave)
  ] : [];

  do_conv_study ? [
    (0.6, model_squares[i], wave) for i=1:length(model_squares)
  ] : [];
]
function mkcase(n, model, bcs)
  (;n, model, bcs)
end
cases = [
  mkcase(n, model, bcs)
  for (n, model, bcs) in paramlist
]

outs1 = []
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
    fun=out.Pih,
    name=name*"/Pi",
    contourargs=(;levels=[2.0^i for i in -1:1:3]|>x->vcat(-x,x)),
    surfaceargs=(;
      colormap=:turbo
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

if do_study
  begin
    scale = 1.5
    nplot_Nu = 100
    xs_Nu = LinRange(0.0, 1.0, 2*nplot_Nu)[1:end]

    map(enumerate([
      1:3,
      4:6
    ])) do (i_indices, indices)

      map(["left", "bottom"]) do wall


        f = Figure(
          size=scale .* (250, 250),
          figure_padding= (0,10,0,0)
        )
        ax = Axis(
          f[1, 1],
          title=(
            wall * " wall, " * (i_indices==1 ? "uniform heating" : "non-uniform heating")
          ),
          xlabel="position " * (wall=="left" ? "y" : "x"),
          ylabel="Nusselt number Nu",
          xminorticksvisible=true,
          yminorticksvisible=true,
          xticks=LinearTicks(5),
          yticks=LinearTicks(5),
          limits=((0.0, 1.0), (0,20))
        )

        map(outs1[indices]) do o
          nb = get_normal_vector(o.btrian)
          NUU = Interpolable(
            # get (-dT/dr) with a transformation to polar coordinates
            (- ∇(o.Th) ⋅ (
              wall=="left" ?
                x->VectorValue([-1.0,0.0]) :
                x->VectorValue([0.0,1.0])
            )); searchmethod=KDTreeSearch(num_nearest_vertices=500)
          )
          cache = return_cache(NUU, Gridap.Point(0.0, 0.0))
          Nu = [
            evaluate!(cache, NUU,
              wall=="left" ?
                Gridap.Point([0.0, x]) :
                Gridap.Point([x, 0.0])
            )
            for x in xs_Nu
          ]
          @info sum(Nu[1:end-1])/length(Nu[1:end-1])
          lines!(
            xs_Nu,
            Nu,
            label="n = $(o.n)",
          )
          CSV.write("nusselt_number_square_$(wall)_$(i_indices)_$(o.n).csv",
            DataFrame(pos = xs_Nu, Nu = Nu)
          )
        end
        axislegend(
          labelsize=10,
          position=:lt,
          orientation=:horizontal,
          nbanks=1,
        )
        f |> display
        save("nusselt_number_square_$(wall)_$(i_indices).pdf", f)
        # export as CSV
      end
    end
  end
end

# Postprocessing of convergence study
if do_conv_study
  begin
    results = DataFrame(
      id = Int[],
      n_elems = Int[],
      n_nodes = Int[],
      Nu_bot_avg = Float64[],
      Sth_max = Float64[],
      Sfl_max = Float64[]
    )

    n_plot = 100
    xs = LinRange(0, 1, n_plot)
    ys = LinRange(0, 1, n_plot)

    map(enumerate([(1:5) .+ (do_study ? 6 : 0)])) do (i_indices, indices)

      map(enumerate(outs1[indices])) do (i, o)

        # mesh stats
        n_elems = length(o.model.grid.cell_type)
        @show n_elems
        n_nodes = length(o.model.grid_topology.vertex_coordinates)  # w. Dirichlet
        @show n_nodes

        # avg(Nu)
        nb = get_normal_vector(o.btrian)
        NUU = Interpolable(
          (- ∇(o.Th) ⋅ (x->VectorValue([0.0,1.0])));
            searchmethod=KDTreeSearch(num_nearest_vertices=500)
        )
        cache = return_cache(NUU, Gridap.Point(0.0, 0.0))
        Nu_inner = [
          evaluate!(cache, NUU, Gridap.Point([x, 0.0]))
          for x in xs
        ]
        Nu_bot_avg = sum(Nu_inner[1:end-1])/length(Nu_inner[1:end-1])
        @show Nu_bot_avg

        # max(Sth)
        search_method = KDTreeSearch(num_nearest_vertices=500)
        Sth_int = Interpolable(o.Sth; searchmethod=search_method)
        cache = return_cache(Sth_int, Gridap.Point(0.0, 0.0))
        function fun_help1(x)
          return evaluate!(cache, Sth_int, Gridap.Point(x))
        end
        zs = [fun_help1([x, y]) for x in xs, y in ys]
        Sth_max = extrema(zs)[2]
        @show Sth_max

        # max(Sfl)
        search_method = KDTreeSearch(num_nearest_vertices=500)
        Sfl_int = Interpolable(o.Sfl; searchmethod=search_method)
        cache = return_cache(Sfl_int, Gridap.Point(0.0, 0.0))
        function fun_help2(x)
          return evaluate!(cache, Sfl_int, Gridap.Point(x))
        end
        zs = [fun_help2([x, y]) for x in xs, y in ys]
        Sfl_max = extrema(zs)[2]
        @show Sfl_max

        push!(
          results,
          (
              id         = i,
              n_elems    = n_elems,
              n_nodes    = n_nodes,
              Nu_bot_avg = Nu_bot_avg,
              Sth_max    = Sth_max,
              Sfl_max    = Sfl_max
          )
      )
      end
    end
    # Display the final table
    @show results

    # Export to CSV
    CSV.write("conv_study_square.csv", results)
  end
end



# 2)
model_annulus = GmshDiscreteModel(
  haskey(ENV, "GITHUB_ACTIONS") ?
  joinpath("../meshes/2.5/co-annulus_unstructured_anisotrop_1.msh") :
  joinpath("../meshes/2.5/co-annulus_unstructured_anisotrop_5.msh")
)
model_annuluses = [GmshDiscreteModel(
  haskey(ENV, "GITHUB_ACTIONS") ?
  joinpath("../meshes/2.5/co-annulus_unstructured_anisotrop_1.msh") :
  joinpath("../meshes/2.5/co-annulus_unstructured_anisotrop_$(i).msh")
) for i in 1:5]
uniform_annulus = (;
  T_diri_tags=["inner", "outer"],
  T_diri_expressions=[1.0, 0.0],
  T_natural_tags=["heatfunction_zero_points"],
  V_diri_tags=["all"]
)
wave_annulus = (;
  T_diri_tags=["inner", "outer"],
  T_diri_expressions=[x->0.5*(sin(atan(x[2],x[1]))+1), 0.0],
  T_natural_tags=["heatfunction_zero_points"],
  V_diri_tags=["all"]
)

paramlist2 = [
  # n, model, bcs, Sfl_lvls
  do_study ? [
    (0.6, model_annulus, uniform_annulus, [10.0^i for i=-2:1:2])
    (1.0, model_annulus, uniform_annulus, [4.0^i for i in -0:1:3])
    (1.4, model_annulus, uniform_annulus, [4.0^i for i in -0:1:3])
    (0.6, model_annulus, wave_annulus, [10.0^i for i=-2:1:2])
    (1.0, model_annulus, wave_annulus, [4.0^i for i in -2:1:3])
    (1.4, model_annulus, wave_annulus, [4.0^i for i in -2:1:3])
  ] : [];
  do_conv_study ? [
      (0.6, model_annuluses[i], wave_annulus, [10.0^i for i=-2:1:2])
      for i=1:length(model_annuluses)
  ] : [];
]
function mkcase(n, model, bcs, Sfl_lvls)
  (;n, model, bcs, Sfl_lvls)
end
cases = [
  mkcase(n, model, bcs, Sfl_lvls)
  for (n, model, bcs, Sfl_lvls) in paramlist2
]

outs2 = []
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
    fun=out.Pih,
    name=name*"/Pi",
    contourargs=(;levels=[2.0^i for i in -1:1:3]|>x->vcat(-x,x)),
    surfaceargs=(;colormap=:turbo)
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

# Nusselt number plots
if do_study
  begin
    scale = 1.5
    nplot_Nu = 100
    phis_Nu = LinRange(0, 2*pi, 2*nplot_Nu)[1:end]

    map(enumerate([1:3,4:6])) do (i_indices, indices)
      map(["inner", "outer"]) do wall

        f = Figure(
          size=scale .* (250, 250),
          figure_padding= (0,11,0,0)
        )
        ax = Axis(
          f[1, 1],
          title = (
            wall * " wall, " * (i_indices==1 ? "uniform" : "non-uniform") *
            " heating"
          ),
          xlabel="angle ϕ",
          ylabel="Nusselt number Nu",
          xminorticksvisible=true,
          yminorticksvisible=true,
          xticks = (0:1/2*pi:2*pi, ["0", "π/2", "π", "3π/2", "2π"]),
          yticks=LinearTicks(5),
          limits=((0.0, 2*pi), (-0.5,13.5))
        )

        map(outs2[indices]) do o
          nb = get_normal_vector(o.btrian)
          Nu = Interpolable(
            # get (-dT/dr) with a transformation to polar coordinates
            (- ∇(o.Th) ⋅ (x->VectorValue([x[1],x[2]] / sqrt(x[1]^2 + x[2]^2))));
              searchmethod=KDTreeSearch(num_nearest_vertices=50)
          )
          cache = return_cache(Nu, Gridap.Point(0.0, 0.0))
          dist = 0
          ri = 2/3
          ro = 5/3-0.001
          Nu = [
            evaluate!(
              cache, Nu,
              wall=="inner" ?
                Gridap.Point([(ri+dist)*cos(p), (ri+dist)*sin(p)]) :
                Gridap.Point([(ro-dist)*cos(p), (ro-dist)*sin(p)])
              )
            for p in phis_Nu
          ]
          lines!(
            phis_Nu,
            Nu,
            label="n = $(o.n)",
          )
          CSV.write("nusselt_number_annulus_$(wall)_$(i_indices)_$(o.n).csv",
            DataFrame(pos = xs_Nu, Nu = Nu)
          )
        end
        axislegend(
          labelsize=10,
          position=:lt,
          orientation=:horizontal,
          nbanks=1,
        )
        f |> display
        save("nusselt_number_annulus_$(wall)_$(i_indices).pdf", f)
      end
    end
  end
end

# Postprocessing of convergence study
if do_conv_study
  begin
    results = DataFrame(
      id = Int[],
      n_elems = Int[],
      n_nodes = Int[],
      Nu_inner_avg = Float64[],
      Sth_max = Float64[],
      Sfl_max = Float64[]
    )

    n_plot = 100
    # xs = LinRange(0, 1, n_plot)
    # ys = LinRange(0, 1, n_plot)
    dist = 0.0
    ri = 2/3
    ro = 5/3
    epss = 0.01
    rs = LinRange(ri+epss, ro-epss, n_plot)
    phis = LinRange(0, 2*pi, 2*n_plot-1)

    map(enumerate([1:5])) do (i_indices, indices)

      map(enumerate(outs2[indices])) do (i, o)

        # mesh stats
        n_elems = length(o.model.grid.cell_types)
        @show n_elems
        n_nodes = length(o.model.grid_topology.vertex_coordinates)  # w. Dirichlet
        @show n_nodes

        # avg(Nu)
        nb = get_normal_vector(o.btrian)
        Nu = Interpolable(
          # get (-dT/dr) with a transformation to polar coordinates
          (- ∇(o.Th) ⋅ (x->VectorValue([x[1],x[2]] / sqrt(x[1]^2 + x[2]^2))));
            searchmethod=KDTreeSearch(num_nearest_vertices=50)
        )
        cache = return_cache(Nu, Gridap.Point(1.0, 1.0))
        Nu_inner = [
          evaluate!(cache, Nu, Gridap.Point([(ri+dist)*cos(p), (ri+dist)*sin(p)]))
          for p in phis
        ]
        Nu_inner_avg = sum(Nu_inner[1:end-1])/length(Nu_inner[1:end-1])
        @show Nu_inner_avg

        # max(Sth)
        search_method = KDTreeSearch(num_nearest_vertices=500)
        Sth_int = Interpolable(o.Sth; searchmethod=search_method)
        cache = return_cache(Sth_int, Gridap.Point(1.0, 1.0))
        function fun_help1(x)
          return evaluate!(cache, Sth_int, Gridap.Point(x))
        end
        xs = [r * cos(phi) for r in rs, phi in phis]
        ys = [r * sin(phi) for r in rs, phi in phis]
        zs = fun_help1.(broadcast((x, y) -> (x, y), xs, ys))
        Sth_max = extrema(zs)[2]
        @show Sth_max

        # max(Sfl)
        search_method = KDTreeSearch(num_nearest_vertices=500)
        Sfl_int = Interpolable(o.Sfl; searchmethod=search_method)
        cache = return_cache(Sfl_int, Gridap.Point(1.0, 1.0))
        function fun_help2(x)
          return evaluate!(cache, Sfl_int, Gridap.Point(x))
        end
        xs = [r * cos(phi) for r in rs, phi in phis]
        ys = [r * sin(phi) for r in rs, phi in phis]
        zs = fun_help2.(broadcast((x, y) -> (x, y), xs, ys))
        Sfl_max = extrema(zs)[2]
        # Sfl_max = sum(zs)/length(zs)
        @show Sfl_max

        push!(
          results,
          (
              id         = i,
              n_elems    = n_elems,
              n_nodes    = n_nodes,
              Nu_inner_avg = Nu_inner_avg,
              Sth_max    = Sth_max,
              Sfl_max    = Sfl_max
          )
      )
      end
    end
    # Display the final table
    @show results

    # Export to CSV
    CSV.write("conv_study_annulus.csv", results)
  end
end
