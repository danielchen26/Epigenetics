using Flux, DiffEqFlux, DifferentialEquations, Plots

## Setup ODE to optimize
function lotka_volterra(du,u,p,t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end
u0 = [1.0,1.0]
tspan = (0.0,10.0)
p = [1.5,1.0,3.0,1.0]
prob = ODEProblem(lotka_volterra,u0,tspan,p)

# # Verify ODE solution
# sol = solve(prob,Tsit5())
# plot(sol)

# Generate data from the ODE
# sol = solve(prob,Tsit5(),saveat=0.1)
# A = sol[1,:] # length 101 vector
# t = 0:0.1:10.0
# scatter!(t,A)

# Build a neural network that sets the cost as the difference from the
# generated data and 1

p = param([2.2, 1.0, 2.0, 0.4]) # Initial Parameter Vector
function predict_rd() # Our 1-layer neural network
  diffeq_rd(p,prob,Tsit5(),saveat=0.1)[1,:]
end
loss_rd() = sum(abs2,x-2 for x in predict_rd()) # loss function

# Optimize the parameters so the ODE's solution stays near 1

data = Iterators.repeated((), 100)
opt = ADAM(0.1)
cb = function () #callback function to observe training
  display(loss_rd())
  # using `remake` to re-create our `prob` with current parameters `p`
  display(plot(solve(remake(prob,p=Flux.data(p)),Tsit5(),saveat=0.1),ylim=(0,6)))
end
# Display the ODE with the initial parameter values.
cb()
Flux.train!(loss_rd, [p], data, opt, cb = cb)























# model-zoo/other/diffeq/neural_ode.jl
u0 = Float32[2.; 0.]
datasize = 30
tspan = (0.0f0,1.5f0)

function trueODEfunc(du,u,p,t)
    true_A = [-0.1 2.0; -2.0 -0.1]
    du .= ((u.^3)'true_A)'
end

test_f = @ode_def begin
    da = m*a - 3n*b + k1
    db = n*b^2  -k2*a
end m n b k1 k2
p = [1., 1., 2., 4., 5.]
u0 = Float32[2.; 0.]
prob = ODEProblem(test_f,u0,tspan,p)
t = range(tspan[1],tspan[2],length=datasize)
prob = ODEProblem(trueODEfunc,u0,tspan)
ode_data = Array(solve(prob,Tsit5(),saveat=t))

dudt = Chain(x -> x.^3,
             Dense(2,10,tanh),
             Dense(10,2))
ps = Flux.params(dudt)
n_ode = x->neural_ode(dudt,x,tspan,Tsit5(),saveat=t,reltol=1e-7,abstol=1e-9)

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



dudt = Chain(x -> x.^2,
             Dense(3,10,tanh),
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
