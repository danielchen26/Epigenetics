

#  -------------------- SDE -------------------
using DiffEqBiological, DifferentialEquations
using Plots;gr()

rn1d = @reaction_network begin
    (a, ep),           2O + Di ↔ Da
     b,                Da → Da + O
    (gama,th),         Di ↔ Dm
     d,                O → ∅
 end a ep b gama th d

params = [1, 0.2, 1., 0.25, 0.25, 0.22]
@add_constraints rn1d begin
    Di + Da + Dm = 1
end


u0=[100.,100.,0.,0.]
tspan = (0,1e4)
# solve Chemical Langevin SDEs
sprob = SDEProblem(rn1d, u0, tspan, params)
ssol  = solve(sprob, EM(), dt=.0001)
plot(ssol)

# solve JumpProblem using Gillespie's Direct Method
u0=[4,0,0,10]
tspan = (0.,1e6)
dprob = DiscreteProblem(rn1d, u0, tspan, params)
jprob = JumpProblem(dprob, Direct(), rn1d)
jsol = solve(jprob, SSAStepper())
plot(jsol, lw = 0.2)
