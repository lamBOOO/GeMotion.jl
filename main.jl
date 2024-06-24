using Gridap
using Gridap.CellData
using Gridap.Arrays
using CairoMakie
using MakieTeX
using LineSearches: BackTracking, StrongWolfe

using Random
Random.seed!(1234)

function simulate2(
  ;Pr=1.0, Ra=1.0, n=1.0, n_elems=20, linesearch=BackTracking(), levels=(;psi=5)
  , T_diri_tags=["leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint", "toprightpoint"], T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0], bcname="turan"
)
  @info Pr, Ra, n, T_diri_tags, T_diri_expressions

  # n_elems = 100
  domain = (0,1,0,1)
  partition = (n_elems, n_elems)
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

  if !isempty(T_diri_tags)
    Θ = TestFESpace(model,reffe_T,conformity=:H1,dirichlet_tags=T_diri_tags)
  else
    Θ = TestFESpace(model,reffe_T,conformity=:H1)
  end

  u_noslip = VectorValue(0,0)  # noslip

  U = TrialFESpace(V,u_noslip)
  P = TrialFESpace(Q)
  T = TrialFESpace(Θ,T_diri_expressions)

  Y = MultiFieldFESpace([V, Q, Θ])
  X = MultiFieldFESpace([U, P, T])

  degree = order
  Ωₕ = Triangulation(model)
  dΩ = Measure(Ωₕ,degree)

  g = VectorValue([0,1])

  γ = 1E-3 # influences convergence of Newton and avoids singularities
  γdot(∇u) = γ + 2*(D(∇u) ⊙ D(∇u))
  D(∇u) = 1/2 * (∇u + ∇u')
  μ(∇u,n) = (x->x^((n-1)/2)) ∘ γdot(∇u)
  dμ(∇u,d∇u,n) = ((n-1)/2) * ((x->x^((n-3)/2)) ∘ γdot(∇u)) * 4 *(D(∇u) ⊙ D(d∇u))

  conv(u,∇u) = (∇u')⋅u
  dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

  a((u,p,T),(v,q,θ)) = ∫(
    0 * Pr*∇(v)⊙∇(u)  # TODO: Make if statement for n==1?
    - (∇⋅v)*p + q*(∇⋅u) + ∇(T) ⋅ ∇(θ) - Ra*Pr*(T)*(g⋅v)
  )dΩ

  b(u,v) = ∫( 2*Pr*μ(∇(u),n)*D(∇(v))⊙D(∇(u)) )dΩ
  db(u,du,v) = ∫(
    2*Pr*(dμ(∇(u),∇(du),n)*D(∇(u)) + μ(∇(u),n)*D(∇(du))) ⊙ D(∇(v))
  )dΩ

  c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
  dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

  d(u,T,θ) = ∫( (u ⋅ ∇(T))*θ )dΩ
  dd(u,du,T,dT,θ) = ∫( ((du ⋅ ∇(T))*θ + (u ⋅ ∇(dT))*θ) )dΩ

  res((u,p,T),(v,q,θ)) = (
    a((u,p,T),(v,q,θ))
    + b(u,v)
    + c(u,v) + d(u,T,θ)
  )
  jac((u,p,T),(du,dp,dT),(v,q,θ)) = (
    a((du,dp,dT),(v,q,θ))
    + db(u,du,v)
    + dc(u,du,v) + dd(u,du,T,dT,θ)
  )

  op = FEOperator(res,jac,X,Y)

  nls = NLSolver(
    show_trace=true, method=:newton
    # , beta=0.1
    , linesearch=linesearch
    # , linesearch=StrongWolfe()
    , ftol=1E-8
    , xtol=1E-10
    )
  solver = FESolver(nls)

  # random initial guess
  # xu = rand(Float64,num_free_dofs(U))
  # # uh0 = FEFunction(U,xu)
  # xp = rand(Float64,num_free_dofs(P))
  # # ph0 = FEFunction(P,xp)
  # xT = rand(Float64,num_free_dofs(T))
  # # Th0 = FEFunction(T,xT)
  # init_guess = FEFunction(X, vcat(xu,xp,xT))
  # result = solve!(init_guess,solver,op)
  # uh = result[1][1]
  # ph = result[1][2]
  # Th = result[1][3]

  # zero initial guess
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

  writevtk(Ωₕ,"results_$(Ra)_$(Pr)_$(n)_$(bcname).vtu",cellfields=[
    "uh"=>uh,"ph"=>ph,"Th"=>Th, "psih"=>psih
  ])

  # plotting
  n_plot = 100
  xs = LinRange(0, 1, n_plot)
  ys = LinRange(0, 1, n_plot)
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
  @info findmax(zs)

  contour!(xs, ys, zs
    ,levels=levels.psi
    ,labels=true
    ,labelsize = 15
    ,color=:black
  )
  save("./streamfunction_$(Ra)_$(Pr)_$(n)_$(bcname).pdf", f)

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
  save("./temperature_$(Ra)_$(Pr)_$(n)_$(bcname).pdf", f)

  # Nusselt number
  # TODO: Add left wall
  nb = get_normal_vector(btrian)
  Nu = Interpolable(
    (∇(Th)⋅nb); searchmethod = KDTreeSearch(num_nearest_vertices=5)
  )

  # entropies
  Sth = interpolate(
    ∇(Th) |> x->inner(x,x)
  , TestFESpace(model,ReferenceFE(lagrangian,Float64,2),conformity=:H1)
  )
  Sfl = interpolate(∇(uh) |> x->1E-4 *(
    2*(inner(x .* x,TensorValue(1,0,0,1)))
    +
    ((inner(x,TensorValue(0,1,1,0))) |> x->x*x)
  ), TestFESpace(model,ReferenceFE(lagrangian,Float64,2),conformity=:H1))
  writevtk(Ωₕ,"entropy_$(Ra)_$(Pr)_$(n)_$(bcname).vtu",cellfields=[
    "Sth"=>Sth, "Sfl" => Sfl
  ])

  xs = LinRange(0, 1, n_plot)
  ys = LinRange(0, 1, n_plot)

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
    # ,levels=10
    ,labels=true
    ,labelsize = 15
    ,color=:black
  )
  save("./entropy_heat_$(Ra)_$(Pr)_$(n)_$(bcname).pdf", f)

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
    # ,levels=10
    ,labels=true
    ,labelsize = 15
    ,labelformatter=formatter
    ,color=:black
  )
  save("./entropy_fluid_$(Ra)_$(Pr)_$(n)_$(bcname).pdf", f)

  return (uh=uh, ph=ph, Th=Th, psih=psih, Nu=Nu, Sth=Sth, Sfl=Sfl, btrian=btrian, model=model, Ωₕ=Ωₕ, Pr=Pr, Ra=Ra)
end

# Define boundary conditions
wave(x) = sin(pi*x[1])
basak_uniform = (;bcname="basak_uniform", T_diri_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"], T_diri_expressions=[1.0,0.0,0.0,0.5,0.5])
basak_wave = (;bcname="basak_wave", T_diri_tags=["botline", "leftline", "rightline", "botleftpoint", "botrightpoint"], T_diri_expressions=[wave,0.0,0.0,0.0,0.0])
turan = (;bcname="turan", T_diri_tags=["leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint", "toprightpoint"], T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0])

# out = simulate2(Pr=0.7, Ra=1E3, n=1.0, 50, linesearch=BackTracking(), levels=(;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))); turan...)

cases=
[
  [0.7, 1E3, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_uniform],
  [0.7, 5*1E3, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.15,0.5,1,1.3] |> x->vcat(x,-x))), basak_uniform],
  [0.7, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,13] |> x->vcat(x,-x))), basak_uniform],
  [0.1, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,4,7,9] |> x->vcat(x,-x))), basak_uniform],
  [1.0, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_uniform],
  [10.0, 1E5, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_uniform],
  [0.015, 1E3, 1.0, 100, BackTracking(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_uniform],
  [0.7, 1E3, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_wave],
  [0.7, 5*1E3, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.15,0.5,1,1.3] |> x->vcat(x,-x))), basak_wave],
  [0.7, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,13] |> x->vcat(x,-x))), basak_wave],
  [0.1, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,4,7,9] |> x->vcat(x,-x))), basak_wave],
  [1.0, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_wave],
  [10.0, 1E5, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([1,5,10,14] |> x->vcat(x,-x))), basak_wave],
  [0.015, 1E3, 1.0, 100, StrongWolfe(), (;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))), basak_wave],
  # # turan
  # [1E3, 1E4, 0.6, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([i for i=1:2:13] |> x->vcat(x,-x))), turan],
  # [1E3, 1E4, 1.0, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:4.5],[5.0]) |> x->vcat(x,-x))), turan],
  # [1E3, 1E4, 1.8, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([i for i=0.1:0.2:1.1] |> x->vcat(x,-x))), turan],
  # [1E3, 1E5, 0.6, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([4,12,20,26,28,29] |> x->vcat(x,-x))), turan],
  # [1E3, 1E5, 1.0, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=1:2:11],[10]) |> x->vcat(x,-x))), turan],
  # [1E3, 1E5, 1.8, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([i for i=0.2:0.4:3.0] |> x->vcat(x,-x))), turan],
  # [1E3, 1E6, 0.6, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=([10,30,50,60,65] |> x->vcat(x,-x))), turan],
  # [1E3, 1E6, 1.0, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=2:4:18],[19]) |> x->vcat(x,-x))), turan],
  # [1E3, 1E6, 1.8, 200, BackTracking(), (;T=[0.1*i for i=1:10],psi=(vcat([i for i=0.5:1.0:5.5],[6]) |> x->vcat(x,-x))), turan],
]

outs = []
for case in cases
  out = simulate2(Pr=case[1], Ra=case[2], n=case[3], n_elems=case[4], linesearch=case[5], levels=case[6]; case[7]...)
  push!(outs, out)
end


# Nusselt number
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
1+1
