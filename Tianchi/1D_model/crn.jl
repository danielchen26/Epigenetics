using DiffEqBiological
using Plots;gr()
using DataFrames, Queryverse



rn1d = @reaction_network begin
    (a, ep),           2O + Di ↔ Da
     b,                Da → Da + O
    (gama,th),         Di ↔ Dm
     d,                O → ∅
 end a ep b gama th d

p = [1., 0.1, 1., 1., 1., 0.5]
@add_constraints rn1d begin
    Di + Da + Dm = 1
end
ss = steady_states(rn1d,p)
stability(ss,rn1d,p)

bif_grid_dia = bifurcation_grid_diagram(rn1d, p, :ep, 0.:0.01:0.25,  :th, (0.1,8.))
plot(bif_grid_dia)



# ----------------------- ODE ------------------------
using DifferentialEquations
module func
    include("../functions.jl")
end

u0 = rand(4);tspan=(1,100.)
prob = ODEProblem(rn1d,u0,tspan,p)
sol = solve(prob,Rosenbrock23())

plot(sol)



# ----- Interact plot for SSS control --------
using Interact
@manipulate for a = 0:0.01:1.0 , ep = 0:0.01:0.5, b = 0:0.01:1.0, th = 0:0.01:0.5, gama = 0:0.01:0.5, d = 0:0.01:1.0
    p = [a,ep,b,gama,th,d]
    ss = steady_states(rn1d,p)
    plot(sort([i[1] for i in ss]), marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
end
