using DifferentialEquations#, DiffEqBiological
using ParameterizedFunctions
using Plots;pyplot(dpi =250)
using Latexify

module epi
    include("functions.jl")
    include("ODE_model_set.jl")
    include("Init.jl")
end




# 1. ============ Single CRN test =============
# randomize Num initial condition for same conservation class
Num = 200
DO_tot = epi.randfixsum(Num,4)
DT_tot = epi.randfixsum(Num,4)
# ------ CRN initial condition ---------
@show u0, p = epi.crn_init(DO_tot,DT_tot,2)
tspan = (0.0,1.5e3)
# ------ Control(choose one) -------
# 1. Oct4 overexpress
ts, cb = epi.make_cb([500,7e2],13,0.06,15,0.36)
# 2. TET overexpress
ts, cb = epi.make_cb([500,7e2],13,0.06,14,0.36,15,0.06)
# 3. Nanog overexpress
ts, cb = epi.make_cb([500,7e2],13,0.36,15,0.06)
# ------ Solving
prob = ODEProblem(epi.paper_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),callback=cb, tstops=ts)
plot(sol, vars = [:N, :T, :O],  title = "CRN Model", ylabel= "Concentration", lw = 1.5, ylims = (0,1.5))







# 2. ============ Randomized initial condition CRN test ============
OTN_SS = []； OTN_SS_all = []
anim = @animate for i  = 1:Num
    # ------- initial condition -----
    u0, p = epi.crn_init(DO_tot,DT_tot,i)
    tspan = (0.0,1e4)

    # ------ Control(choose one) -------
    # 1. Oct4 overexpress
    # ts, cb = epi.make_cb([5e3,7e3],13,0.06,15,0.36)
    # 2. TET overexpress
    ts, cb = epi.make_cb([5e3,7e3],13,0.06,14,0.36,15,0.06)
    # 3. Nanog overexpress
    # ts, cb = epi.make_cb([5e3,7e3],13,0.36,15,0.06)

    # ------ Solving -------
    prob = ODEProblem(epi.paper_ode,u0,tspan,p)
    sol = solve(prob, DynamicSS(Rodas5()),callback=cb, tstops=ts)
    push!(OTN_SS, sol[end][[1,2,9]]); push!(OTN_SS_all, sol[end])
    # ------ Plotting ------
    plot(sol, vars = [:N, :T, :O],  w = 1.5, ylims = (0,1.5),
    title = "CRN with random initial $i", xlabel = "Time", ylabel= "Concentration")
end
gif(anim, "/Users/chentianchi/Desktop/Oct4_random_init.gif", fps = 5)







# ===================================   test reduced model
paper_reduce_ode = @ode_def begin
    dN = -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1
    dT = (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    dO = (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3


p = [0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6]
tspan = (0.0,1e2)
prob_rd = ODEProblem(paper_reduce_ode,[0.,0.1,0.],tspan,p)
sol_rd = solve(prob_rd, DynamicSS(Rodas5()))
plot(sol_rd, w = 1.5, ylims = (0,1.), title = "Reduced Model N T O", xlabel = "Time", ylabel= "Concentration")

for i in collect(1:10)
    prob_rd = ODEProblem(paper_reduce_ode,rand(3),tspan,p)
    sol_rd = solve(prob_rd, DynamicSS(Rodas5()))
    plot!(sol_rd)
end




display(latexify(paper_reduce_ode.symjac))




# # ================ parameters estimation =================
# using DiffEqParamEstim, DiffEqBayes
# using RecursiveArrayTools, LeastSquaresOptim
# import Optim
# function f2(du,u,p,t)
#   du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
#   du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
# end
#
# u0 = [1.0;1.0]
# tspan = (0.0,10.0)
# p = [1.5,1.0,3.0,1.0]
# prob = ODEProblem(f2,u0,tspan,p)
# sol = solve(prob,Tsit5())
# t = collect(range(0,stop=10,length=200))
#  # for VectorOfArray
# randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
# data = convert(Array,randomized)
#
# cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),
#                                       maxiters=10000,verbose=false)
# result_bfgs = Optim.optimize(cost_function, [1.3,0.8,2.8,1.2], Optim.BFGS())
#
# cost_function = build_lsoptim_objective(prob,t,data,Tsit5())
# x = [1.3,0.8,2.8,1.2]
# res = optimize!(LeastSquaresProblem(x = x, f! = cost_function,
#                 output_length = length(t)*length(prob.u0)),
#                 LeastSquaresOptim.Dogleg(),LeastSquaresOptim.LSMR())
