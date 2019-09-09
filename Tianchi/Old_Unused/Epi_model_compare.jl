using DifferentialEquations, DiffEqBiological
using Latexify, DataFrames
using Plots

# ----------- Importing Models -------------
#################################################

module Model2
include("Animation_gen.jl")
include("Model2.jl")
include("Model2_params_db.jl")
end

module Model3
include("Animation_gen.jl")
include("Model3.jl")
include("Model3_params_db.jl")
end




# -------------  Show imported Models parameters -----------
################################################################

show(Model2.params, allrows = true)
show(Model3.params, allrows = true)




# -------------- Model initialization ---------------------
################################################################

# ---- Model2
p2 = Model2.params.values
# ---- Model3
p3 = Model3.params.values

# ------ Make sure parameters in the Model and in the DataFrame are in the same order
show(Model2.params.names)         # DataFrame parameters
show(Model2.O_N_auto.params)      # Model parameters


tspan = (0., 1e3)
u0 = zeros(Float64,15)
u0[[3,6,10,13]] .= 0.3       # Check below which variables should turn on

# Unmethylated promoters turned on initially
# ------ Model2
display(Model2.O_N_auto.syms_to_ints)
show(Model2.O_N_auto.syms[[3,4,6,10,13,14]])
show(Model2.O_N_auto.syms)
display(Model2.O_N_auto.params_to_ints)
# ------ Model3
display(Model3.O_N_auto.syms_to_ints)
show(Model3.O_N_auto.syms[[3,6,10,13,14]])
show(Model3.O_N_auto.params)


# Show ODE model from CRN
show(Model2.O_N_auto.f_symfuncs)





# ------------------ solve ODEs ---------------------
#######################################################

#  ----------  Model2  ---------
oprob = ODEProblem(Model2.O_N_auto, u0, tspan, p2)
osol  = solve(oprob, Rodas5())
plot(osol, vars = getindex(Model2.O_N_auto.syms, [1,6,11]), legend = true)
#  ----------  Model3  ---------
oprob = ODEProblem(Model3.O_N_auto, u0, tspan, p3)
osol  = solve(oprob, Tsit5())
plot(osol, vars = getindex(Model3.O_N_auto.syms, [1,6,11]), legend = true)




#  ---------------- Make Animation --------------------
############################################################

gr(dpi = 200)

# 1. Simulation with [:D, :Dₒ, :TET, :D♇, :Dₙ, :Dᴺ] turned on
# ------- Model2 Animation
anim = Model2.anim_gen(u0,p2,tspan,Model2.O_N_auto,1,10,[3,4,6,10,13,14],[1], 100, [1,6,11])
gif(anim, "/Users/chentianchi/Desktop/demethy_strength_variation2.gif", fps = 2)
# ------- Model3 Animation
anim = Model3.anim_gen(u0,p3,tspan,Model3.O_N_auto,0.1,1,[3,4,6,10,13,14],[1,6,11],1,1)
gif(anim, "/tmp/demethy_strength_variation3.gif", fps = 2)




# 2. Simulation with [:D̄, :TET] turned on
# ------- Model2 Animation
anim = Model2.anim_gen(u0,p2,tspan,Model2.O_N_auto,0.5,10,[5,6],[1,6,11],1,1)
gif(anim, "/tmp/Initial_methylated.gif", fps = 2)
# ------- Model3 Animation
anim = Model3.anim_gen(u0,p3,tspan,Model3.O_N_auto,0.1,1,[5,6],[1,6,11],1,1)
gif(anim, "/tmp/Initial_methylated.gif", fps = 5)









# # ---------------- Model2 ODE model (CL changed) ------------- testing consistent with CRN
# Model2.O_N_auto_ode.syms
# # Model2.O_N_auto_ode
# p = Model2.params.values
# u0 = zeros(Float64,15)
# u0[[3,6,10,13]] .= 0.3
# tspan = (0.0,1e3)
# prob = ODEProblem(Model2.O_N_auto_ode,u0,tspan,p)
# sol = solve(prob)
# # plotly()
# plot(sol, vars = [1,6,11])
#
# getindex(Model2.O_N_auto_ode.syms, [1,6,11])
#
#
# pyplot()
# anim = @animate for i=0:1:10
#     u0[[3,4,6,10,13,14]] .= i
#     u0[1] = 100
#     # getindex(Model2.O_N_auto_ode.syms, [3,4,6,10,13,14])
#     # explicit ode
#     prob = ODEProblem(Model2.O_N_auto_ode,u0,tspan,p)
#     sol = solve(prob, Rosenbrock23())
#     plot(sol, vars = getindex(Model2.O_N_auto_ode.syms, [1,6,11]), legend = true, lw =1)
#     title!("Model Evolution : Strength Response : DM_initial = $i ")
#     xaxis!("Time")
#     yaxis!("Concentration")
# end
# gif(anim, "/Users/chentianchi/Desktop/ode_test_M2.gif", fps = 2)
# # ---------------- Model2 ODE model ------------- testing
































##  --------##  -------##  --------##  -------##  --------##  -------##
##  -----------------    Some testing    ----------------------------##
##  --------##  -------##  --------##  -------##  --------##  -------##



# Make Animation
# pyplot()
anim = @animate for i=0:0.2:10
    u0[[3,4,6,10,13,14]] .= i
    oprob = ODEProblem(Model2.O_N_auto, u0, tspan, p2)
    # steady states solutions
    sssol = solve(oprob,DynamicSS(Tsit5()))
    plot(sssol, vars = getindex(Model2.O_N_auto.syms, [1,6,11]),  legend = true, line=(1.2,:dash)) # tspan=(0.0,1.0),

    # # regular solution
    # osol  = solve(oprob, Tsit5())
    # plot(osol, vars = getindex(O_N_auto.syms, [1,6,11]), legend = true, line=(2,:dash))
    ylims!(0,1)
    title!("Model Evolution : Strength Response : DM_initial = $i ")
    xaxis!("Time")
    yaxis!("Concentration")
end
gif(anim, "/tmp/demethylation_strength_variation.gif", fps = 2)












# ---------------   Testing the ode control example
using DifferentialEquations, ParameterizedFunctions,Plots;gr()
mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::Array{T,1}
end

ts = [1e2,3e2]
condition(u,t,integrator) = t in ts
function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p[25] = 1.0
  elseif integrator.t == ts[2]
      integrator.p[25] = 0.0
  end
end
cb = DiscreteCallback(condition, affect!, save_positions=(true,true));

epi_model = @ode_def_bare begin # m1 -> O     m2 -> N     m3 -> TET
    dO = -O*β - O*δₒ - O^2*a₀ + 2*O₂*d₀ + k*Dₒ*α₁ + m1
    dO₂ = Dᴺ*ζ₀ + (1/2)*O^2*a₀ - O₂*d₀ + rᵦ*D♇ + r₀*Dₒ - β*O₂ - D*r₁*O₂ - O₂*Dₙ*ζ₁ - r₂*O₂*Dₜ
    dD = -D*aₐ + r₀*Dₒ + β*Dₕ - D*r₁*O₂
    dDₒ = -r₀*Dₒ + D*r₁*O₂
    dD̄ = D*aₐ + d₂*C₂ - NT*D̄*k₃ - TET*D̄*a₂
    dTET = D♇*ϕ₁ + NT*η₀ - TET*δₜ + d₂*C₂ + k₂*C₂ - β*TET - N*TET*η₁ - TET*D̄*a₂ + m3
    dC₂ = -d₂*C₂ - k₂*C₂ - β*C₂ + TET*D̄*a₂
    dDₕ = k₂*C₂ - β*Dₕ + NT*D̄*k₃
    dDₜ = rᵦ*D♇ - NT*k₄*Dₜ - r₂*O₂*Dₜ
    dD♇ = -rᵦ*D♇ + NT*k₄*Dₜ + r₂*O₂*Dₜ
    dN = -N*β - N*δₒ - N^2*a₀ + NT*η₀ + 2*N₂*d₀ - N*TET*η₁ + k*Dᴺ*α₁ + m2
    dN₂ = (1/2)*N^2*a₀ - N₂*d₀ + r₀*Dᴺ - β*N₂ - r₁*N₂*Dₙ
    dDₙ =  Dᴺ*ζ₀ + r₀*Dᴺ - O₂*Dₙ*ζ₁ - r₁*N₂*Dₙ
    dDᴺ = -Dᴺ*ζ₀ - r₀*Dᴺ + O₂*Dₙ*ζ₁ + r₁*N₂*Dₙ
    dNT = -NT*η₀ + N*TET*η₁
end a₀ d₀ r₁ r₀ α₀ α₁ k β δₒ aₐ a₂ d₂ k₂ ϕ₀ ϕ₁ r₂ rᵦ δₜ ζ₁ ζ₀ η₁ η₀ k₃ k₄ m1 m2 m3




# Solving
u00 = zeros(Float64,15); u00[[3,4,6,10,13,14]] .= 0.5
p00 = Array([p3;[0,0,0]])
u0 = SType(u00, p00)
tspan = (0.0,1e3)
p = p00
prob = ODEProblem(epi_model,u0,tspan,p)

sol = solve(prob, DynamicSS(Rosenbrock23()),callback=cb, tstops=ts)
plot(sol, vars = getindex(epi_model.syms, [1,6,11]),  legend = true, line=(1.5,:dash))













#  Animation
u00 = zeros(Float64,15); u00[[3,4,6,10,13,14]] .= 0.5
tspan = (0.0,1e3)
anim = @animate for i=0:1e3:1e5
#  control
    ts = [1e2,3e2]
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
      if integrator.t == ts[1]
          integrator.p[25] = 1.0
      elseif integrator.t == ts[2]
          integrator.p[25] = 0.0
      end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));

#   solve
    p3[24] = i
    p00 = Array([p3;[0,0,0]])
    u0 = SType(u00, p00)
    p = p00
    prob = ODEProblem(epi_model,u0,tspan,p)
    # steady states solutions
    sssol = solve(prob,DynamicSS(Rodas5()), callback=cb,tstops=ts)
    plot(sssol, vars = getindex(epi_model.syms, [1,6,11]),  legend = true, line=(1.2,:dash))
    # ylims!(0,1)
    title!("Model Evolution : Strength Response : DM_initial = $i ")
    xaxis!("Time")
    yaxis!("Concentration")
end

gif(anim, "/Users/chentianchi/Desktop/test.gif", fps = 2)














































# # working example-----------------------------------------------
# ts = [5,10]
# condition(u,t,integrator) = t in ts
# function affect!(integrator)
#   if integrator.t == ts[1]
#       integrator.p[1] = 10.0
#   elseif integrator.t == ts[2]
#       integrator.p[1] = 0.0
#   end
# end
#
# save_positions = (true,true)
# cb = DiscreteCallback(condition, affect!, save_positions=save_positions);
#
# test = @ode_def_bare test_model begin
#     dx = -0.5x + m1
#     dv = -0.5v + m2
# end m1 m2
#
# u0 = SType([10.0;10.0], [0.0,1.0])
# tspan = (0.0,15.0)
# p = [0.0,1.0]
# prob = ODEProblem(test,u0,tspan,p)
# sol = solve(prob,Tsit5(),callback=cb,tstops=ts)
# plot(sol)
