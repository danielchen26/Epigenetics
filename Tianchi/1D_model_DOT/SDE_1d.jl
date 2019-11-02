#  -------------------- SDE -------------------
using DiffEqBiological, DifferentialEquations
using Plots;gr()
using LinearAlgebra
rn1d = @reaction_network begin
    (a, ep),           2O + Di ↔ Da
     b,                Da → Da + O
     be,               Di → Di + O
    (gama,th),         Di ↔ Dm
     d,                O → ∅
 end a ep b be gama th d



rn1d_leak = @reaction_network begin
    (α, ϵ),            2O + Di ↔ Da
     β,                Da → Da + O
     η,                Di → Di + O
    (γ,θ),             Di ↔ Dm
     δ,                O → ∅
 end α ϵ β η γ θ δ
# params = [1., 0.2, 1., 0.25, 0.25, 0.22]
# params = [.5, 0.25, 0.5, 0.49, 0.05, 7.39] # 0 ~ 4
# params = [5., 10., 5.06, 10., 0.04, 6.69] # 0 ~ 58
# params = [10., 3.05, 10., 10., 0.04, 1.75]
params = [5., 5., 5., 0.03, 5., 5., 5.]
@add_constraints rn1d_leak begin
    Di + Da + Dm = 5.
end
ss = steady_states(rn1d_leak,params)
# stability(ss,rn1d,params)

# ---- below is the correct stability -------
sort!(ss, by = x -> x[1])
@show ss
function CRN_Eig(rn,params,ss)
    J=rand(length(rn.syms),length(rn.syms))
    t=0
    jacfun(rn)(J,ss,params,t)
    return eigvals(J), eigvecs(J)
end
@show [round.(CRN_Eig(rn1d,params,i)[1],digits = 2) for i in ss]


# ======== Interact =============
using Interact
@manipulate for a = 0:0.01:10.0 , ep = 0:0.01:10., b = 0:0.01:10.0, eta = 0:0.01:3.0, gama = 0:0.01:10.,  th= 0:0.01:10., d = 0:0.01:10.0
    p = [a,ep,b,eta,gama,th,d]
    ss = steady_states(rn1d_leak,p)
    plot(sort([i[1] for i in ss]), marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
    # println("The SS states are:")
    # sort!(ss, by = x -> x[1])
    # @show ss
    # @show [round.(EigenVector(rn1d,p,i),digits = 2) for i in ss]
    #
    # eigs = [JE_stability(i,rn1d,p)[2] for i in ss]
    # println("The eigs are:")
    # # @show eigs
    # @show hcat(eigs...)
end





# ------- Try this params ------
params = [.5, 0.25, 0.5, 0.49, 0.05, 7.39] # 0 ~ 4
params = [5., 10., 5.06, 10., 0.04, 6.69] # 0 ~ 58
# --------- ODE -------------
@manipulate for O = 0:200., Di = 0:100., Dm = 0:100., tmax = 0.:0.1:1e2
    u0 = [O,Di,Dm,100. - Di - Dm];tspan=(1e-9,tmax)
    prob = ODEProblem(rn1d,u0,tspan,params)
    sol = solve(prob,Rosenbrock23())
    plot(sol)
end

# solve Chemical Langevin SDEs

u0=[1e2,1e3,0.,0.]
tspan = (0,1e2)
sprob = SDEProblem(rn1d, u0, tspan, params)
ssol  = solve(sprob, SimpleTauLeaping(), dt=.1 )
plot(ssol)

# solve JumpProblem using Gillespie's Direct Method
u0 = [10.,1.,3.,1.]
tspan = (0.,1e3)
params = [5., 5., 5., 5., 0.03, 5., 5.]
dprob = DiscreteProblem(rn1d_leak, u0, tspan, params)
jprob = JumpProblem(dprob, Direct(), rn1d_leak)
jsol = solve(jprob, SSAStepper())
plot(jsol,vars = [:O], lw = 0.2)



using StatsPlots,DataFrames
using DataVoyager
plotly()
df = DataFrame(jsol')
Voyager(jsol)
@df df scatter(:x1, :x2, :x3)




#  ==== Interact SSA ========
@manipulate for a = 0:0.01:1.0 , ep = 0:0.01:0.5, b = 0:0.01:1.0, th = 0:0.01:0.5, gama = 0:0.01:0.5, d = 0:0.01:1.0
    p = [a,ep,b,gama,th,d]
    u0 = [11.,500.,100.,400.]
    tspan = (0.,1e2)
    dprob = DiscreteProblem(rn1d, u0, tspan, p)
    jprob = JumpProblem(dprob, Direct(), rn1d)
    jsol = solve(jprob, SSAStepper())
    plot(jsol) # vars=[:O]
end
