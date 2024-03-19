
using Gridap
n = 100
domain = (0,1,0,1)
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

# labelling:
# 1-4: points in corners
# 5-6-7-8 = bot-top-left-right

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"top",[6,])
add_tag_from_tags!(labels,"5",[5,])
add_tag_from_tags!(labels,"6",[6,])
add_tag_from_tags!(labels,"7",[7,])
add_tag_from_tags!(labels,"8",[8,])
add_tag_from_tags!(labels,"not-top",[1,2,3,4,5,7,8])

D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,labels=labels,dirichlet_tags=["not-top","top"])

reffe_T = ReferenceFE(lagrangian,Float64,order)

reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
Q = TestFESpace(model,reffeₚ,conformity=:L2,constraint=:zeromean)

Θ = TestFESpace(model,reffe_T,conformity=:H1,dirichlet_tags=["5", "6", "7", "8"])

uD0 = VectorValue(0,0)
uD1 = VectorValue(1,0)
U = TrialFESpace(V,[uD0,uD1])
P = TrialFESpace(Q)
T = TrialFESpace(Θ,[1.0,2.0,3.0,4.0])

Y = MultiFieldFESpace([V, Q, Θ])
X = MultiFieldFESpace([U, P, T])

degree = order
Ωₕ = Triangulation(model)
dΩ = Measure(Ωₕ,degree)

const Re = 10.0
conv(u,∇u) = Re*(∇u')⋅u
dconv(du,∇du,u,∇u) = conv(u,∇du)+conv(du,∇u)

a((u,p,T),(v,q,θ)) = ∫( ∇(v)⊙∇(u) - (∇⋅v)*p + q*(∇⋅u) + ∇(T) ⋅ ∇(θ) )dΩ

c(u,v) = ∫( v⊙(conv∘(u,∇(u))) )dΩ
dc(u,du,v) = ∫( v⊙(dconv∘(du,∇(du),u,∇(u))) )dΩ

res((u,p,T),(v,q,θ)) = a((u,p,T),(v,q,θ)) + c(u,v)
jac((u,p,T),(du,dp,dT),(v,q,θ)) = a((du,dp,dT),(v,q,θ)) + dc(u,du,v)

op = FEOperator(res,jac,X,Y)

using LineSearches: BackTracking
nls = NLSolver(
  show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)

uh, ph, Th = solve(solver,op)

writevtk(Ωₕ,"ins-results",cellfields=["uh"=>uh,"ph"=>ph,"Th"=>Th])
