using Gridap
using Gridap.CellData
using Gridap.Arrays
using CairoMakie
using MakieTeX
using LineSearches: BackTracking, StrongWolfe

function simulate2(;Pr=1.0, Ra=1.0, bc=:one, linesearch=BackTracking(), levels=(;psi=5))
  @info Pr, Ra, bc

  n = 40
  domain = (0,1,0,1)
  partition = (n,n)
  model = CartesianDiscreteModel(domain,partition)
  btrian = BoundaryTriangulation(model)

  # labelling:
  # 1-2-3-4 = botleftpoint-botrightpoint-topleftpoint-toprightpoint
  # 5-6-7-8 = botline-topline-leftline-rightline

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"botleftpoint",[1,])
  add_tag_from_tags!(labels,"botrightpoint",[2,])
  add_tag_from_tags!(labels,"topleftpoint",[3,])
  add_tag_from_tags!(labels,"toprightpoint",[4,])
  add_tag_from_tags!(labels,"botline",[5,])
  add_tag_from_tags!(labels,"topline",[6,])
  add_tag_from_tags!(labels,"leftline",[7,])
  add_tag_from_tags!(labels,"rightline",[8,])
  add_tag_from_tags!(labels,"all",[1,2,3,4,5,6,7,8])

  order = 2
  reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["all"])

  reffe_T = ReferenceFE(lagrangian,Float64,order)

  reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
  Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)

  Θ = TestFESpace(model,reffe_T,conformity=:H1,dirichlet_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"])

  u_noslip = VectorValue(0,0)  # noslip
  wave(x) = sin(pi*x[1])

  if bc == :one
    bc_fun = 1.0
  elseif  bc == :wave
    bc_fun = wave
  else
    error("bc not supported")
  end

  U = TrialFESpace(V,u_noslip)
  P = TrialFESpace(Q)
  T = TrialFESpace(Θ,[bc_fun,0.0,0.0,0.5,0.5])

  Y = MultiFieldFESpace([V, Q, Θ])
  X = MultiFieldFESpace([U, P, T])

  degree = order
  Ωₕ = Triangulation(model)
  dΩ = Measure(Ωₕ,degree)

  # Ra = 1E5
  # Pr = 0.7
  g = VectorValue([0,1])
  conv(u,∇u) = (∇u')⋅u
  dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

  a((u,p,T),(v,q,θ)) = ∫( Pr*∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) + ∇(T) ⋅ ∇(θ) - Ra*Pr*(T)*(g⋅v) )dΩ

  c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
  dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

  d(u,T,θ) = ∫( (u ⋅ ∇(T))*θ )dΩ
  dd(u,du,T,dT,θ) = ∫( ((du ⋅ ∇(T))*θ + (u ⋅ ∇(dT))*θ) )dΩ

  res((u,p,T),(v,q,θ)) = a((u,p,T),(v,q,θ)) + c(u,v) + d(u,T,θ)
  jac((u,p,T),(du,dp,dT),(v,q,θ)) = a((du,dp,dT),(v,q,θ)) + dc(u,du,v) + dd(u,du,T,dT,θ)

  op = FEOperator(res,jac,X,Y)

  nls = NLSolver(
    show_trace=true, method=:newton
    # , beta=0.1
    , linesearch=linesearch
    # , linesearch=StrongWolfe()
    , ftol=1E-15
    , xtol=1E-10
    )
  solver = FESolver(nls)
  uh, ph, Th = solve(solver,op)

  # stream function
  reffe_psi = ReferenceFE(lagrangian,Float64,1)
  test_psi = TestFESpace(model,reffe_psi,conformity=:H1,dirichlet_tags=["all"])
  trial_psi = TrialFESpace(test_psi,0.0)
  aa(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
  bb(v) = ∫( v * (∇(uh) ⊙ TensorValue(0,-1,1,0)) )*dΩ
  op = AffineFEOperator(aa,bb,trial_psi,test_psi)
  ls = LUSolver()
  solver = LinearFESolver(ls)
  psih = solve(solver,op)

  writevtk(Ωₕ,"results_$(Ra)_$(Pr)_$(bc).vtu",cellfields=["uh"=>uh,"ph"=>ph,"Th"=>Th, "psih"=>psih])

  # plotting
  xs = LinRange(0, 1, n)
  ys = LinRange(0, 1, n)
  search_method = KDTreeSearch(num_nearest_vertices=5)

  f = Figure(
    size = (500, 500)
    ,figure_padding = (0, 10, 0, 0)
  )
  Axis(
    f[1, 1]
    ,aspect=1
    ,title="stream function ψ"
    ,xlabel="x"
    ,ylabel="y"
    ,xminorticksvisible = true
    ,yminorticksvisible = true
    ,xticks = LinearTicks(6)
    ,yticks = LinearTicks(6)
    ,limits = (0, 1, 0, 1)
  )
  psihi = Interpolable(psih; searchmethod=search_method)
  cache_psi = return_cache(psihi, Gridap.Point(0.0, 0.0))
  function helper_psi(x)
    return evaluate!(cache_psi, psihi, Gridap.Point(x))
  end
  zs = [helper_psi([x,y]) for x in xs, y in ys]

  contour!(xs, ys, zs
    ,levels=levels.psi
    ,labels=true
    ,labelsize = 15
    ,color=:black
  )
  save("./streamfunction_$(Ra)_$(Pr)_$(bc).pdf", f)

  Thi = Interpolable(Th; searchmethod=search_method)
  cache_T = return_cache(Thi, Gridap.Point(0.0, 0.0))
  function helper_T(x)
    return evaluate!(cache_T, Thi, Gridap.Point(x))
  end
  zs = [helper_T([x,y]) for x in xs, y in ys]
  f = Figure(
    size = (500, 500)
    ,figure_padding = (0, 10, 0, 0)
  )
  Axis(
    f[1, 1]
    ,aspect=1
    ,title="temperature T"
    ,xlabel="x"
    ,ylabel="y"
    ,xminorticksvisible = true
    ,yminorticksvisible = true
    ,xticks = LinearTicks(6)
    ,yticks = LinearTicks(6)
    ,limits = (0, 1, 0, 1)
  )
  contour!(xs, ys, zs
    ,levels=levels.T
    ,labels=true
    ,labelsize = 15
    ,color=:black
  )
  save("./temperature_$(Ra)_$(Pr)_$(bc).pdf", f)

  # Nusselt number
  # TODO: Add left wall
  nb = get_normal_vector(btrian)
  Nu = Interpolable((∇(Th) ⋅ nb); searchmethod = KDTreeSearch(num_nearest_vertices=5))

  # entropies
  Sth = interpolate(
    ∇(Th) |> x->inner(x,x)
  , TestFESpace(model,ReferenceFE(lagrangian,Float64,2),conformity=:H1))
  Sfl = interpolate(∇(uh) |> x->1E-4 *(
    2*(inner(x .* x,TensorValue(1,0,0,1)))
    +
    ((inner(x,TensorValue(0,1,1,0))) |> x->x*x)
  ), TestFESpace(model,ReferenceFE(lagrangian,Float64,2),conformity=:H1))
  writevtk(Ωₕ,"entropy_$(Ra)_$(Pr)_$(bc).vtu",cellfields=["Sth"=>Sth, "Sfl" => Sfl])

  n = 100
  xs = LinRange(0, 1, n)
  ys = LinRange(0, 1, n)

  cache_Sth = return_cache(Sth, Gridap.Point(0.0, 0.0))
  cache_Sfl = return_cache(Sfl, Gridap.Point(0.0, 0.0))
  z_Sfl = [evaluate!(cache_Sfl, Sfl, Gridap.Point([x,y])) for x in xs, y in ys]
  z_Sth = [evaluate!(cache_Sth, Sth, Gridap.Point([x,y])) for x in xs, y in ys]

  f = Figure(
    size = (500, 500)
    ,figure_padding = (0, 10, 0, 0)
  )
  Axis(
    f[1, 1]
    ,aspect=1
    ,title="local heat entropy S_θ"
    ,xlabel="x"
    ,ylabel="y"
    ,xminorticksvisible = true
    ,yminorticksvisible = true
    ,xticks = LinearTicks(6)
    ,yticks = LinearTicks(6)
    ,limits = (0, 1, 0, 1)
  )
  contour!(xs, ys, z_Sth
    ,levels=vcat(1:1:7, 10:5:30, 0.01,0.1,0.5,0.25)
    ,labels=true
    ,labelsize = 15
    ,color=:black
  )
  save("./entropy_heat_$(Ra)_$(Pr)_$(bc).pdf", f)

  f = Figure(
    size = (500, 500)
    ,figure_padding = (0, 10, 0, 0)
  )
  Axis(
    f[1, 1]
    ,aspect=1
    ,title="local fluid entropy S_ψ"
    ,xlabel="x"
    ,ylabel="y"
    ,xminorticksvisible = true
    ,yminorticksvisible = true
    ,xticks = LinearTicks(6)
    ,yticks = LinearTicks(6)
    ,limits = (0, 1, 0, 1)
  )
  function formatter(level::Real)::String
      lev_short = round(level; digits = 3)
      string(isinteger(lev_short) ? round(Int, lev_short) : lev_short)
  end
  contour!(xs, ys, z_Sfl
    ,levels=vcat(0.005:0.005:0.05)
    ,labels=true
    ,labelsize = 15
    ,labelformatter=formatter
    ,color=:black
  )
  save("./entropy_fluid_$(Ra)_$(Pr)_$(bc).pdf", f)

  return (uh=uh, ph=ph, Th=Th, psih=psih, Nu=Nu, Sth=Sth, Sfl=Sfl, btrian=btrian, model=model, Ωₕ=Ωₕ, Pr=Pr, Ra=Ra)
end

out = simulate2(Pr=10, Ra=1E5, bc=:wave, linesearch=StrongWolfe(), levels=(;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))))

cases=
[
  [0.7, 1E3, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x)))],
  [0.7, 5*1E3, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.15,0.5,1,1.3] |> x->vcat(x,-x)))],
  [0.7, 1E5, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,13] |> x->vcat(x,-x)))],
  [0.1, 1E5, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,4,7,9] |> x->vcat(x,-x)))],
  [1.0, 1E5, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x)))],
  [10.0, 1E5, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x)))],
  [0.015, 1E3, :one, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x)))],
  [0.7, 1E3, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x)))],
  [0.7, 5*1E3, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.15,0.5,1,1.3] |> x->vcat(x,-x)))],
  [0.7, 1E5, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,13] |> x->vcat(x,-x)))],
  [0.1, 1E5, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,4,7,9] |> x->vcat(x,-x)))],
  [1.0, 1E5, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x)))],
  [10.0, 1E5, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x)))],
  [0.015, 1E3, :wave, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x)))],
]

outs = []
for case in cases
  out = simulate2(Pr=case[1], Ra=case[2], bc=case[3], linesearch=case[4], levels=case[5])
  push!(outs, out)
end


# Nusselt number
begin
begin
nn = 100
xs = 0 |> x->LinRange(0+x, 1-x, nn)
# update_theme!(
#   palette = (color=Makie.wong_colors(), marker = [:circle, :xcross, :pentagon, :diamond],),
#   ScatterLines = (cycle = Cycle([:color, :marker], covary = true),)
# )
update_theme!(
  palette = (color=Makie.wong_colors(), marker = [:circle, :xcross, :pentagon, :diamond],),
  Lines = (cycle = Cycle([:color, :marker], covary = true),)
)
f = Figure(
  # size = (500, 500)
  figure_padding = (0, 10, 0, 0)
)
ax = Axis(
  f[1, 1]
  ,aspect=300/180
  ,title="local Nusselt number Nu, bottom wall"
  ,xlabel="x"
  ,ylabel="y"
  ,xminorticksvisible = true
  ,yminorticksvisible = true
  ,xticks = LinearTicks(5)
  ,yticks = LinearTicks(3)
  ,limits = ((0.1, 0.9), (0,15))
)

map(outs[[1,3,6]]) do o
  cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
  lines!(xs, [evaluate!(cache, o.Nu, Gridap.Point([x,0])) for x in xs]
    ,ticks = LinearTicks(8)
    # ,color = :black
    ,label = "Pr=$(o.Pr), Ra=$(o.Ra)"
  )
end
map(outs[[8,10,13]]) do o
  cache = return_cache(o.Nu, Gridap.Point(0.0, 0.0))
  lines!(xs, [evaluate!(cache, o.Nu, Gridap.Point([x,0])) for x in xs]
    ,ticks = LinearTicks(8)
    # ,color = :black
    ,label = "Pr=$(o.Pr), Ra=$(o.Ra)"
    ,linestyle=:dash
  )
end
axislegend(labelsize=10)
save("./nusselt_bot.pdf", f)
f
end
end
1+1
