using Gridap
using Gridap.CellData
using Gridap.Arrays
using CairoMakie
using LineSearches: BackTracking, StrongWolfe

using Random

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
  name="sim",
  Pr=1.0,
  Ra=1.0,
  n=1.0,
  model=CartesianDiscreteModel((0, 1, 0, 1), (10, 10)),
  nlsolver_opts=(;
    show_trace=true,
    method=:newton,
    linesearch=BackTracking(),
    ftol=1E-8,
    xtol=11E-6
  ),
  nlsolver_custom_init_guess=Float64[],
  nlsolver_init_guess_type=:zero,
  jac_scaling=1.0,
  T_diri_tags=[
    1,2,3,4,5,6
  ],
  T_diri_expressions=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
  V_diri_tags=[1, 2, 3, 4, 5, 6, 7, 8],
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
    model, reffeᵤ, conformity=:H1, labels=labels, dirichlet_tags=V_diri_tags
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
    + b(u, v)  # disable for newtonian
    + c(u, v) + d(u, T, θ)
  )
  jac((u, p, T), (du, dp, dT), (v, q, θ)) = jac_scaling * (
    a((du, dp, dT), (v, q, θ))
    + db(u, du, v)  # disable for newtonian
    + dc(u, du, v) + dd(u, du, T, dT, θ)
  )

  "Setup operator" |> println
  op = FEOperator(res, jac, X, Y)
  # use AD, see https://gridap.github.io/Tutorials/stable/pages/t004_p_laplacian/
  # op = FEOperator(res, X, Y)

  println("Setup nonlinear solver")
  nls = NLSolver(; nlsolver_opts...)
  solver = FESolver(nls)

  println("Solve nonlinear problem")
  if nlsolver_init_guess_type == :zero
    println("zero initial guess")
    result = solve(solver, op)
  elseif nlsolver_init_guess_type == :ones
    println("constant one initial guess")
    xu = ones(Float64,num_free_dofs(U))
    xp = ones(Float64,num_free_dofs(P))
    xT = ones(Float64,num_free_dofs(T))
    # uh0 = FEFunction(U,xu)
    # ph0 = FEFunction(P,xp)
    # Th0 = FEFunction(T,xT)
    init_guess = FEFunction(X, vcat(xu,xp,xT))
    result = solve!(init_guess,solver,op)[1]
  elseif nlsolver_init_guess_type == :random
    println("random initial guess")
    Random.seed!(1234)
    xu = rand(Float64, num_free_dofs(U))
    xp = rand(Float64, num_free_dofs(P))
    xT = rand(Float64, num_free_dofs(T))
    init_guess = FEFunction(X, vcat(xu, xp, xT))
    result = solve!(init_guess, solver, op)[1]
  elseif nlsolver_init_guess_type == :custom
    println("custom initial guess")
    init_guess = FEFunction(X, nlsolver_custom_init_guess)
    result = solve!(init_guess,solver,op)[1]
  end
  uh, ph, Th = result

  # stream function
  reffe_psi = ReferenceFE(lagrangian, Float64, 1)
  test_psi = TestFESpace(
    model, reffe_psi, conformity=:H1, dirichlet_tags=V_diri_tags
  )
  trial_psi = TrialFESpace(test_psi, 0.0)
  aa(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
  bb(v) = ∫(v * (∇(uh) ⊙ TensorValue(0, -1, 1, 0))) * dΩ
  op = AffineFEOperator(aa, bb, trial_psi, test_psi)
  ls = LUSolver()
  solver = LinearFESolver(ls)
  psih = solve(solver, op)

  fname = joinpath(name, "results.vtu")
  println("Write results to $(fname)")
  writevtk(Ωₕ, fname, cellfields=[
    "uh" => uh, "ph" => ph, "Th" => Th, "psih" => psih
  ])

  return (
    uh=uh,
    ph=ph,
    Th=Th,
    result=result,
    psih=psih,
    btrian=btrian,
    model=model,
    Ωₕ=Ωₕ,
    Pr=Pr,
    Ra=Ra
  )
end
simulate = simulate3
