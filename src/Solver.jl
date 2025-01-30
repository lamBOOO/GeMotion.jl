using Gridap
using Gridap.CellData
using Gridap.Arrays
using CairoMakie
using MakieTeX
using LineSearches: BackTracking, StrongWolfe

using Random
Random.seed!(1234)

export simulate

"""
    simulate(;name="sim",Pr=1.0,Ra=1.0,n=1.0,model=CartesianDiscreteModel((0, 1, 0, 1), (20, 20)),nlsolver_opts=(;show_trace=true,method=:newton,linesearch=BackTracking(),ftol=1E-8,xtol=1E-10),T_diri_tags=["leftline","rightline","botleftpoint","botrightpoint","topleftpoint","toprightpoint"],T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0])

Simulate the flow in a square cavity with a moving lid and a temperature gradient.

# Arguments
- `name::String`: Name of the simulation.
- `Pr::Float64`: Prandtl number.
- `Ra::Float64`: Rayleigh number.
- `n::Float64`: Power of the viscosity.
- `model::DiscreteModel`: Number of elements in each direction.
- `nlsolver_opts::Dict`: Options for the nonlinear solver.
- `T_diri_tags::Array{String}`: Tags for the Dirichlet boundary conditions for the temperature.
- `T_diri_expressions::Array{Float64}`: Expressions for the Dirichlet boundary conditions for the temperature.

# Returns
- `uh::FEFunction`: Velocity field.
- `ph::FEFunction`: Pressure field.
- `Th::FEFunction`: Temperature field.
- `psih::FEFunction`: Stream function.
- `Nu::Interpolable`: Nusselt number.
- `Sth::FEFunction`: Local heat entropy.
- `Sfl::FEFunction`: Local fluid entropy.
- `btrian::BoundaryTriangulation`: Boundary triangulation.
- `model::CartesianDiscreteModel`: Model.
- `Ωₕ::Triangulation`: Triangulation.
- `Pr::Float64`: Prandtl number.
- `Ra::Float64`: Rayleigh number.

# Example
```julia
simulate(;name="sim",Pr=1.0,Ra=1.0,n=1.0,model=CartesianDiscreteModel((0, 1, 0, 1), (20, 20)),nlsolver_opts=(;show_trace=true,method=:newton,linesearch=BackTracking(),ftol=1E-8,xtol=1E-10),T_diri_tags=["leftline","rightline","botleftpoint","botrightpoint","topleftpoint","toprightpoint"],T_diri_expressions=[0.0,1.0,0.0,1.0,0.0,1.0])
```

"""
function simulate3(;
  name = "sim",
  Pr = 1.0,
  Ra = 1.0,
  n = 1.0,
  model = CartesianDiscreteModel((0, 1, 0, 1), (20, 20)),
  nlsolver_opts = (;
    show_trace = true,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1E-8,
    xtol = 11E-6
  ),
  T_diri_tags = [
    "leftline", "rightline", "botleftpoint", "botrightpoint", "topleftpoint",
    "toprightpoint"
  ],
  T_diri_expressions = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
)
  @info Pr, Ra, n, T_diri_tags, T_diri_expressions
  @info "Start"

  @info "Output folder: $(joinpath(pwd(),name))"
  mkpath(name)

  btrian = BoundaryTriangulation(model)
  labels = get_face_labeling(model)

  order = 2
  reffeᵤ = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
  V = TestFESpace(
    model, reffeᵤ, conformity=:H1, labels=labels, dirichlet_tags=["all"]
  )

  reffe_T = ReferenceFE(lagrangian, Float64, order)

  reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space=:P)
  Q = TestFESpace(model, reffeₚ, conformity=:L2, constraint=:zeromean)

  if !isempty(T_diri_tags)
    Θ = TestFESpace(model, reffe_T, conformity=:H1, dirichlet_tags=T_diri_tags)
  else
    Θ = TestFESpace(model, reffe_T, conformity=:H1)
  end

  u_noslip = VectorValue(0, 0)  # noslip

  U = TrialFESpace(V, u_noslip)
  P = TrialFESpace(Q)
  T = TrialFESpace(Θ, T_diri_expressions)

  Y = MultiFieldFESpace([V, Q, Θ])
  X = MultiFieldFESpace([U, P, T])

  degree = order
  Ωₕ = Triangulation(model)
  dΩ = Measure(Ωₕ, degree)

  g = VectorValue([0, 1])

  reg = 1E-3 # influences convergence of Newton and avoids singularities
  gamma_dot(∇u) = reg + 2 * (D(∇u) ⊙ D(∇u))
  D(∇u) = 1 / 2 * (∇u + ∇u')
  μ(∇u, n) = (x -> x^((n - 1) / 2)) ∘ gamma_dot(∇u)
  dμ(∇u, d∇u, n) = (
    ((n - 1) / 2) * ((x -> x^((n - 3) / 2)) ∘ gamma_dot(∇u))
    * 4 * (D(∇u) ⊙ D(d∇u))
  )

  conv(u, ∇u) = (∇u') ⋅ u
  dconv(du, ∇du, u, ∇u) = conv(u, ∇du) + conv(du, ∇u)

  a((u, p, T), (v, q, θ)) = ∫(
    0 * Pr * ∇(v) ⊙ ∇(u)  # TODO: Make if statement for n==1?
    -
    (∇ ⋅ v) * p + q * (∇ ⋅ u) + ∇(T) ⋅ ∇(θ) - Ra * Pr * (T) * (g ⋅ v)
  )dΩ

  b(u, v) = ∫(2 * Pr * μ(∇(u), n) * D(∇(v)) ⊙ D(∇(u)))dΩ
  db(u, du, v) = ∫(
    2 * Pr * (dμ(∇(u), ∇(du), n) * D(∇(u)) + μ(∇(u), n) * D(∇(du))) ⊙ D(∇(v))
  )dΩ

  c(u, v) = ∫(v ⊙ (conv ∘ (u, ∇(u))))dΩ
  dc(u, du, v) = ∫(v ⊙ (dconv ∘ (du, ∇(du), u, ∇(u))))dΩ

  d(u, T, θ) = ∫((u ⋅ ∇(T)) * θ)dΩ
  dd(u, du, T, dT, θ) = ∫(((du ⋅ ∇(T)) * θ + (u ⋅ ∇(dT)) * θ))dΩ

  res((u, p, T), (v, q, θ)) = (
    a((u, p, T), (v, q, θ))
    + b(u, v)
    + c(u, v) + d(u, T, θ)
  )
  jac((u, p, T), (du, dp, dT), (v, q, θ)) = (
    a((du, dp, dT), (v, q, θ))
    + db(u, du, v)
    + dc(u, du, v) + dd(u, du, T, dT, θ)
  )

  "Setup operator" |> println
  op = FEOperator(res, jac, X, Y)

  println("Setup nonlinear solver")
  nls = NLSolver(; nlsolver_opts...)
  solver = FESolver(nls)

  # random initial guess
  # println("Solve nonlinear problem")
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
  println("Solve nonlinear problem")
  uh, ph, Th = solve(solver, op)

  # stream function
  reffe_psi = ReferenceFE(lagrangian, Float64, 1)
  test_psi = TestFESpace(
    model, reffe_psi, conformity=:H1, dirichlet_tags=["all"]
  )
  trial_psi = TrialFESpace(test_psi, 0.0)
  aa(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
  bb(v) = ∫(v * (∇(uh) ⊙ TensorValue(0, -1, 1, 0))) * dΩ
  op = AffineFEOperator(aa, bb, trial_psi, test_psi)
  ls = LUSolver()
  solver = LinearFESolver(ls)
  psih = solve(solver, op)

  writevtk(Ωₕ, joinpath(name, "results.vtu"), cellfields=[
    "uh" => uh, "ph" => ph, "Th" => Th, "psih" => psih
  ])

  return (
    uh=uh,
    ph=ph,
    Th=Th,
    psih=psih,
    btrian=btrian,
    model=model,
    Ωₕ=Ωₕ,
    Pr=Pr,
    Ra=Ra
  )
end
simulate = simulate3
