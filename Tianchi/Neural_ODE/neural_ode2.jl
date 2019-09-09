##### -- ============ My code ============ -- #####
using ParameterizedFunctions,DifferentialEquations
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
# plot(sol_rd, w = 1.5, ylims = (0,1.), title = "Reduced Model N T O", xlabel = "Time", ylabel= "Concentration")

u0 = Float32[0.1; 0.5; 0.1]
datasize = 30
tspan = (0.0f0,1.5f0)

t = range(tspan[1],tspan[2],length=datasize)
prob = ODEProblem(paper_reduce_ode,u0,tspan,p)
ode_data = Array(solve(prob,DynamicSS(Rodas5()),saveat=t))



dudt = Chain(Dense(3,10,tanh),
             Dense(10,3))
ps = Flux.params(dudt)
n_ode = x->neural_ode(dudt,x,tspan,DynamicSS(Rodas5()),saveat=t,reltol=1e-7,abstol=1e-9)

pred = n_ode(u0) # Get the prediction using the correct initial condition
scatter(t,ode_data[1,:],label="data")
scatter!(t,Flux.data(pred[1,:]),label="prediction")

scatter(t,ode_data[2,:],label="data")
scatter!(t,Flux.data(pred[2,:]),label="prediction")

function predict_n_ode()
  n_ode(u0)
end
loss_n_ode() = sum(abs2,ode_data .- predict_n_ode())

data = Iterators.repeated((), 1000)
opt = ADAM(0.1)
cb = function () #callback function to observe training
  display(loss_n_ode())
  # plot current prediction against data
  cur_pred = Flux.data(predict_n_ode())
  pl = scatter(t,ode_data[1,:],label="data")
  scatter!(pl,t,cur_pred[1,:],label="prediction")
  display(plot(pl))
end

# Display the ODE with the initial parameter values.
cb()

Flux.train!(loss_n_ode, ps, data, opt, cb = cb)






















# ------ Parameter Estimation and Bayesian Analysis -------
using DiffEqParamEstim, DifferentialEquations
# using DiffEqBayes

function f(du,u,p,t)
  du[1] = dx = p[1]*u[1] - u[1]*u[2]
  du[2] = dy = -3*u[2] + u[1]*u[2]
end

u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5]
prob = ODEProblem(f,u0,tspan,p)


sol = solve(prob,Tsit5())
t = collect(range(0,stop=10,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),
                                     maxiters=10000,verbose=false)

vals = 0.0:0.1:10.0
using Plots; gr()
plot(vals,[cost_function(i) for i in vals],yscale=:log10,
  xaxis = "Parameter", yaxis = "Cost", title = "1-Parameter Cost Function",
  lw = 3)




using Optim
result = optimize(cost_function, 0.0, 10.0)

lower = [0.0]
upper = [3.0]
result = optimize(cost_function, [1.42], lower, upper, Fminbox{BFGS}())









#  we can use the same tools to estimate multiple parameters simultaneously. Let's use the Lotka-Volterra equation with all parameters free:
function f2(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
end

u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(f2,u0,tspan,p)


sol = solve(prob,Tsit5())
t = collect(range(0,stop=10,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)


cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),
                                      maxiters=10000,verbose=false)
result_bfgs = Optim.optimize(cost_function, [1.3,0.8,2.8,1.2], Optim.BFGS())



cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data,differ_weight=0.3,data_weight=0.7),
                                      maxiters=10000,verbose=false)
@show result_bfgs = Optim.optimize(cost_function, [1.3,0.8,2.8,1.2], Optim.BFGS())


using LeastSquaresOptim
x = [1.3,0.8,2.8,1.2]
res = optimize!(LeastSquaresProblem(x = x, f! = cost_function,
                output_length = length(t)*length(prob.u0)),
                LeastSquaresOptim.Dogleg(),LeastSquaresOptim.LSMR())

using Flux, DiffEqFlux
function delay_lotka_volterra(du,u,h,p,t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = (α - β*y)*h(p,t-0.1)[1]
  du[2] = dy = (δ*x - γ)*y
end
h(p,t) = ones(eltype(p),2)
prob = DDEProblem(delay_lotka_volterra,[1.0,1.0],h,(0.0,10.0),constant_lags=[0.1])

p = param([2.2, 1.0, 2.0, 0.4])
params = Flux.Params([p])
function predict_rd_dde()
  Array(diffeq_rd(p,prob,MethodOfSteps(Tsit5()),saveat=0.1))
end
loss_rd_dde() = sum(abs2,x-1 for x in predict_rd_dde())
loss_rd_dde()

function lotka_volterra_noise(du,u,p,t)
  du[1] = 0.1u[1]
  du[2] = 0.1u[2]
end
prob = SDEProblem(lotka_volterra,lotka_volterra_noise,[1.0,1.0],(0.0,10.0))

p = param([2.2, 1.0, 2.0, 0.4])
params = Flux.Params([p])
function predict_fd_sde()
  diffeq_fd(p,Array,202,prob,SOSRI(),saveat=0.1)
end
loss_fd_sde() = sum(abs2,x-1 for x in predict_fd_sde())
loss_fd_sde()

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cb = function ()
  display(loss_fd_sde())
  display(plot(solve(remake(prob,p=Flux.data(p)),SOSRI(),saveat=0.1),ylim=(0,6)))
end

# Display the ODE with the current parameter values.
cb()

Flux.train!(loss_fd_sde, params, data, opt, cb = cb)


















using DiffEqParamEstim, DifferentialEquations
function f2(du,u,p,t)
  du[1] = dx = p[1]*u[1] - p[2]*u[1]*u[2]
  du[2] = dy = -p[3]*u[2] + p[4]*u[1]*u[2]
end

u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(f2,u0,tspan,p)




sol = solve(prob,Tsit5())
t = collect(range(0,stop=10,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)



using Optim
cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),
                                      maxiters=10000,verbose=false)
@show result_bfgs = Optim.optimize(cost_function, [1.3,0.8,2.8,1.2], Optim.BFGS())


using LeastSquaresOptim
cost_function = build_lsoptim_objective(prob,t,data,Tsit5())

x = [1.3,0.8,2.8,1.2]
res = optimize!(LeastSquaresProblem(x = x, f! = cost_function,
                output_length = length(t)*length(prob.u0)),
                LeastSquaresOptim.Dogleg(),LeastSquaresOptim.LSMR())
println(res.minimizer)






function f(du,u,p,t)
  dx = p[1]*u[1] - u[1]*u[2]
  dy = -3*u[2] + u[1]*u[2]
end

u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5]
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob,Tsit5())

t = collect(range(0,stop=10,length=200))
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

obj = build_loss_objective(prob,Tsit5(),L2Loss(t,data),maxiters=10000)

using NLopt
opt = Opt(:LN_COBYLA, 1)
min_objective!(opt, obj)
(minf,minx,ret) = NLopt.optimize(opt,[1.3])








f1 = function (du,u,p,t)
  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
  du[2] = -3.0 * u[2] + u[1]*u[2]
end
p = [1.5,1.0]
u0 = [1.0;1.0]
tspan = (0.0,10.0)
prob1 = ODEProblem(f1,u0,tspan,p)
sol = solve(prob1,Tsit5())

using RecursiveArrayTools # for VectorOfArray
t = collect(range(0,stop=10,length=200))
function generate_data(sol,t)
  randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
  data = convert(Array,randomized)
end
aggregate_data = convert(Array,VectorOfArray([generate_data(sol,t) for i in 1:100]))
using Distributions
distributions = [fit_mle(Normal,aggregate_data[i,j,:]) for i in 1:2, j in 1:200]
using DiffEqParamEstim
obj = build_loss_objective(prob1,Tsit5(),LogLikeLoss(t,distributions),
                                     maxiters=10000,verbose=false)


using Plots; gr()
rg = collect(range(0.5,stop=5,length=100))
heatmap(rg,rg,[obj([j,i]) for i in rg, j in rg],
     yscale=:log10,xlabel="Parameter 1",ylabel="Parameter 2",
     title="Likelihood Landscape")

plot(rg,[obj([i,1.0]) for i in rg],lw=3,
      title="Parameter 1 Likelihood (Parameter 2 = 1.5)",
      xlabel = "Parameter 1", ylabel = "Objective Function Value")



using BlackBoxOptim
bound1 = Tuple{Float64, Float64}[(0.5, 5),(0.5, 5)]
result = bboptimize(obj;SearchRange = bound1, MaxSteps = 11e3)
