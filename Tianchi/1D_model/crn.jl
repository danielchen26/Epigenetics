using DiffEqBiological
using Plots;gr()
using DataFrames, Queryverse



rn1d = @reaction_network begin
    (a, ep),           2O + Di ↔ Da
     b,                Da → Da + O
    (gama,th),         Di ↔ Dm
     d,                O → ∅
 end a ep b gama th d

p = [1.,1.,1.,1.,1.,1.]
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




using Interact
using Makie
@manipulate for ep = 0:0.01:0.5
    p[2] = ep
    prob = ODEProblem(rn1d,u0,tspan,p)
    sol = solve(prob,Rosenbrock23())
    plot(sol)
end


@manipulate for ep = 0:0.01:0.5
    p[2] = ep
    ss = steady_states(rn1d,p)
    scatter(ss[1][1],ss[1][1])
end


# need 1d vis for SSS
