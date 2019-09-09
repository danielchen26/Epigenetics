using DifferentialEquations
using ParameterizedFunctions

test_f = @ode_def begin
    dx = a*x - 3b*y
    dy = b*y + 2b*x
end a b
p = [1,2]
test_f.jac()
test_f.paramjac()
test_f.symjac


test_f.jac(rand(2),p)
test_f.jac(p, rand(2))



tspan = (0.0,1e2)
prob_rd = ODEProblem(test_f,[0.,0.1],tspan,p)
sol_rd = solve(prob_rd, DynamicSS(Rodas5()))

# Working example
test_f = @ode_def begin
    da = m*a - 3n*b + k1
    db = n*b^2 + v^2*a -k2
end m n b k1 v k2
p = [1,6,2,3,4,5]
t = 0.0
out = rand(2,2)
test_f.jac(out,rand(2),p,t)
@show out
