using DiffEqBiological, LinearAlgebra
using Plots;gr()
using DataFrames, Queryverse, Latexify



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
    # plot(sort([i[1] for i in ss]), marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
    println("The SS states are:")
    @show [i[1] for i in ss]
    eigs = [JE_stability(i,rn1d,p)[2] for i in ss]
    println("The eigs are:")
    # @show eigs
    @show hcat(eigs...)
end

test = ([JE_stability(i,rn1d,p)[2] for i in ss])
hcat(test...)





# -------- Jacobian -----------
function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,p,t)
    return (jac,eigen(jac).values)
end

@show jac, eig = JE_stability(ss[end],rn1d,p)
eig





# ----------- 1d distribution of sss -------------
using VegaLite
using DifferentialEquations

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

u0 = [rand(1); 1.0; 0.0; 0.0];tspan = (1,100.)
prob = ODEProblem(rn1d,u0,tspan,p)
sol = solve(prob,Rosenbrock23())

ss_set = []
x_set = []
for  x = 0:0.01:1.
    push!(x_set,x)

    u0 = [x; 1.0; 0.0; 0.0];tspan = (1,1000.)
    prob = ODEProblem(rn1d,u0,tspan,p)
    sol = solve(prob,Rosenbrock23())
    push!(ss_set,sol[end][1])
end


df = DataFrame(id = x_set,O = ss_set)
df|>
@vlplot(
    :point,
    x=:id,
    y=:O,
    width=400,
    height=400,
    background =:white
)






#  ----heat -----

using Interact
@manipulate for a = 0:0.01:1.0 , ep = 0:0.01:0.5, b = 0:0.01:1.0, th = 0:0.01:0.5, d = 0:0.01:1.0


    p = [1, 0.2, 1., 0.25, 0.25, 0.22]
    Xs = []
    γs = []
    SSs = []
    for γ = 0.: 0.01:0.9, x = 0:0.01:0.9
        push!(Xs,x);push!(γs,γ)
        # p[4] = γ;
        p = [a,ep,b,γ,th,d]
        u0 = [x; 1.0; 0.0; 0.0];tspan = (1,1000.)
        prob = ODEProblem(rn1d,u0,tspan,p)
        sol = solve(prob,Rosenbrock23())
        v =  abs(sol[end][1] - 4) < abs(sol[end][1] - 0) ? 1 : 0
        # @show v
        push!(SSs,v)
    end



    df = DataFrame(X = convert(Array{Float64,1},Xs), γ = convert(Array{Float64,1},γs), st =convert(Array{Float64,1},SSs))
    df|>
    @vlplot(
        :rect,
        x={:X, bin={maxbins=100}},
        y={:γ, bin={maxbins=100}},
        color ="st:o",
        background = :white
    )

end

sol

extrema(sol[end])
sol[end][1]






#  ------ control γ --------
p = [1, 0.2, 1., 0.25, 0.25, 0.22]
Xs = []
γs = []
SSs = []
for γ = 0.: 0.01:0.9, x = 0:0.01:0.9
    push!(Xs,x);push!(γs,γ)
    p[4] = γ;
    u0 = [x; 0.0; 0.0; 1.0];tspan = (1,1000.)
    prob = ODEProblem(rn1d,u0,tspan,p)
    sol = solve(prob,Rosenbrock23())
    @show sol[end]
    v =  abs(sol[end][1] - 4) < abs(sol[end][1] - 0) ? 1 : 0
    # @show v
    push!(SSs,v)
end



df = DataFrame(X = convert(Array{Float64,1},Xs), γ = convert(Array{Float64,1},γs), st =convert(Array{Float64,1},SSs))
df|>
@vlplot(
    :rect,
    x={:X, bin={maxbins=100}},
    y={:γ, bin={maxbins=100}},
    color ="st:o",
    background = :white
)



#  ------ control θ --------
p = [1, 0.2, 1., 0.25, 0.25, 0.22]
Xs = []
θs = []
SSs = []
for θ = 0.: 0.01 :3, x = 0:0.01:3
    push!(Xs,x);push!(θs,θ)
    p[5] = θ;
    u0 = [x; 0.0; 0.0; 1.0];tspan = (1,1000.)
    prob = ODEProblem(rn1d,u0,tspan,p)
    sol = solve(prob,Rosenbrock23())
    @show sol[end]
    v =  abs(sol[end][1] - 4) < abs(sol[end][1] - 0) ? 1 : 0
    # @show v
    push!(SSs,v)
end



df = DataFrame(X = convert(Array{Float64,1},Xs), θ = convert(Array{Float64,1},θs), st =convert(Array{Float64,1},SSs))
df|>
@vlplot(
    :rect,
    x={:X, bin={maxbins=100}},
    y={:θ, bin={maxbins=100}},
    color ="st:o",
    width=300, height=200,
    background = :white
)
