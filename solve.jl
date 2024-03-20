using Gridap
import Random
Random.seed!(1234)

n = 100
domain = (0,1,0,1)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

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
U = TrialFESpace(V,u_noslip)
P = TrialFESpace(Q)
T = TrialFESpace(Θ,[1.0,0.0,0.0,0.5,0.5])

Y = MultiFieldFESpace([V, Q, Θ])
X = MultiFieldFESpace([U, P, T])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

# const Re = 1.0
# const Raa = 1E3
# const Prr = 0.7
Raa = 1E5
Prr = 0.7
const g = VectorValue([0,1])
conv(u,∇u) = (∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

# a((u,p,T),(v,q,θ)) = ∫( Prr*∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) + ∇(T) ⋅ ∇(θ) - Raa*Prr*(g⋅∇(T))*(g⋅v) )dΩ
a((u,p,T),(v,q,θ)) = ∫( Prr*∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) + ∇(T) ⋅ ∇(θ) - Raa*Prr*(T)*(g⋅v) )dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

d(u,T,θ) = ∫( (u ⋅ ∇(T))*θ )dΩ
dd(u,du,T,dT,θ) = ∫( ((du ⋅ ∇(T))*θ + (u ⋅ ∇(dT))*θ) )dΩ

res((u,p,T),(v,q,θ)) = a((u,p,T),(v,q,θ)) + c(u,v) + d(u,T,θ)
jac((u,p,T),(du,dp,dT),(v,q,θ)) = a((du,dp,dT),(v,q,θ)) + dc(u,du,v) + dd(u,du,T,dT,θ)

op = FEOperator(res,jac,X,Y)

using LineSearches: BackTracking


# # initial guess nonzero
# xu = rand(Float64,num_free_dofs(U))
# # uh0 = FEFunction(U,xu)

# xp = rand(Float64,num_free_dofs(P))
# # ph0 = FEFunction(P,xp)

# xT = rand(Float64,num_free_dofs(T))
# # Th0 = FEFunction(T,xT)

# init_guess = FEFunction(X, vcat(xu,xp,xT))

# out = solve!(init_guess,solver,op)

# uh = out[1][1]
# ph = out[1][2]
# Th = out[1][3]

nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking()
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



writevtk(Ωₕ,"ins-results",cellfields=["uh"=>uh,"ph"=>ph,"Th"=>Th, "psih"=>psih])


# plot
using CairoMakie
using MakieTeX
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

xs = LinRange(0, 1, n)
ys = LinRange(0, 1, n)

using Gridap.CellData
using Gridap.Arrays
search_method = KDTreeSearch(num_nearest_vertices=5)

psihi = Interpolable(psih; searchmethod=search_method)
const cache5 = return_cache(psihi, Gridap.Point(0.0, 0.0))
function helper4(x)
  return evaluate!(cache5, psihi, Gridap.Point(x))
end
zs = [helper4([x,y]) for x in xs, y in ys]

# zs = [psih(VectorValue([x,y])) for x in xs, y in ys]

contour!(xs, ys, zs
  ,levels=vcat([1,5,10,13],-[1,5,10,13])
  ,labels=true
  # , transparancy = true
  ,labelsize = 15
  # ,labelfont = :bold
  # ,colormap=:grays
  ,color=:black
  # ,labelcolor=:black
)
save("./streamfunction.pdf", f)
@show f

Thi = Interpolable(Th; searchmethod=search_method)
const cache6 = return_cache(Thi, Gridap.Point(0.0, 0.0))
function helper5(x)
  return evaluate!(cache6, Thi, Gridap.Point(x))
end
zs = [helper5([x,y]) for x in xs, y in ys]

f = Figure(
  size = (500, 500)
  ,figure_padding = (0, 10, 0, 0)
)
Axis(
  f[1, 1]
  ,aspect=1
  ,title="temperature θ"
  ,xlabel="x"
  ,ylabel="y"
  ,xminorticksvisible = true
  ,yminorticksvisible = true
  ,xticks = LinearTicks(6)
  ,yticks = LinearTicks(6)
  ,limits = (0, 1, 0, 1)
)
contour!(xs, ys, zs
  ,levels=[i*0.1 for i=1:10]
  ,labels=true
  # , transparancy = true
  ,labelsize = 15
  # ,labelfont = :bold
  # ,colormap=:grays
  ,color=:black
  # ,labelcolor=:black
)
save("./temperature.pdf", f)
@show f
