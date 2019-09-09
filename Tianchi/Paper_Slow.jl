using DifferentialEquations, DiffEqBiological
using ParameterizedFunctions
using Plots;pyplot(dpi =250)


# ------- constraint dissosiation rate
Kₒ  = 0.3          # Kₒ is the OCT4 dissociation constant
Kₙₜ =  0.2          # KNT is the [N|T] complex dissociation constant
Kd =  0.1          # Kd dimerization constant for NT

ps = @reaction_network begin
    # N T complex
    (a₁, Kd*a₁),         N + T ↔ NT
    (d, d),             NT + NT ↔ NT2

    # --- Promoter binding
    # Oct4
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₀ᴼ ↔ D₁₀ᴼ
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₁ᴼ ↔ D₁₁ᴼ
    (aₒ, Kₒ*aₒ),          O + D₀₀ᴼ ↔ D₀₁ᴼ
    (aₒ, Kₒ*aₒ),          O + D₁₀ᴼ ↔ D₁₁ᴼ
    # TET
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₀ᵀ ↔ D₁₀ᵀ
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₁ᵀ ↔ D₁₁ᵀ
    (aₒ, Kₒ*aₒ),          O + D₀₀ᵀ ↔ D₀₁ᵀ
    (aₒ, Kₒ*aₒ),          O + D₁₀ᵀ ↔ D₁₁ᵀ
    # Nanog
    (aₒ, Kₒ*aₒ),          O + D₀ᴺ ↔ D₁ᴺ

    # --- Protein production
    αₜ,                 D₁₁ᵀ → D₁₁ᵀ + T
    αₒ,                 D₁₁ᴼ → D₁₁ᴼ + O
    αₙ,                 D₁ᴺ → D₁ᴺ + N

    # --- Dilution and Degradation
    (δ,δ,δ),           (N, T, O) → ∅
    #(β,β),             (NT, NT2) → ∅
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β


p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0.]
u0 = zeros(15); u0[2] = 0.5; u0[5] =1; u0[10] =1;u0[14] =1;
tspan = (0.,1e2)

oprob = ODEProblem(ps, u0, tspan, p)
osol  = solve(oprob, Rosenbrock23())
crn0 = plot(osol, vars = [:N], color = [:red],title = "Paper Model", ylabel= "Concentration")
crn0 = plot!(osol, vars = [:T], color = [:green],title = "Paper Model", ylabel= "Concentration")
crn0 = plot!(osol, vars = [:O], color = [:blue],title = "Paper Model", ylabel= "Concentration")

savefig(crn0,"/Users/chentianchi/Desktop/epi_plots/crn0.png")








# -------------- CRN -> ODE model with parameters control
mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::Array{T,1}
end
ts = [5e2,7e2]
condition(u,t,integrator) = t in ts

paper_ode = @ode_def begin
    dN = -N*δ + αₙ*D₁ᴺ - N*T*a₁ + NT*Kd*a₁ + m1
    dT = -T*δ + αₜ*D₁₁ᵀ - N*T*a₁ + NT*Kd*a₁ + m2
    dNT = -NT*β - d*NT^2 + 2*d*NT2 + N*T*a₁ - NT*Kd*a₁
    dNT2 = (1/2)*d*NT^2 - d*NT2 - β*NT2 - NT2*aₙₜ*D₀₀ᴼ - NT2*aₙₜ*D₀₀ᵀ - NT2*aₙₜ*D₀₁ᴼ - NT2*aₙₜ*D₀₁ᵀ + aₙₜ*Kₙₜ*D₁₀ᴼ + aₙₜ*Kₙₜ*D₁₀ᵀ + aₙₜ*Kₙₜ*D₁₁ᴼ + aₙₜ*Kₙₜ*D₁₁ᵀ
    dD₀₀ᴼ = Kₒ*aₒ*D₀₁ᴼ - NT2*aₙₜ*D₀₀ᴼ - O*aₒ*D₀₀ᴼ + aₙₜ*Kₙₜ*D₁₀ᴼ
    dD₁₀ᴼ = Kₒ*aₒ*D₁₁ᴼ + NT2*aₙₜ*D₀₀ᴼ - O*aₒ*D₁₀ᴼ - aₙₜ*Kₙₜ*D₁₀ᴼ
    dD₀₁ᴼ = -Kₒ*aₒ*D₀₁ᴼ - NT2*aₙₜ*D₀₁ᴼ + O*aₒ*D₀₀ᴼ + aₙₜ*Kₙₜ*D₁₁ᴼ
    dD₁₁ᴼ = -Kₒ*aₒ*D₁₁ᴼ + NT2*aₙₜ*D₀₁ᴼ + O*aₒ*D₁₀ᴼ - aₙₜ*Kₙₜ*D₁₁ᴼ
    dO = -O*δ + αₒ*D₁₁ᴼ + Kₒ*aₒ*D₀₁ᴼ + Kₒ*aₒ*D₀₁ᵀ + Kₒ*aₒ*D₁ᴺ + Kₒ*aₒ*D₁₁ᴼ + Kₒ*aₒ*D₁₁ᵀ - O*aₒ*D₀ᴺ - O*aₒ*D₀₀ᴼ - O*aₒ*D₀₀ᵀ - O*aₒ*D₁₀ᴼ - O*aₒ*D₁₀ᵀ + m3
    dD₀₀ᵀ = Kₒ*aₒ*D₀₁ᵀ - NT2*aₙₜ*D₀₀ᵀ - O*aₒ*D₀₀ᵀ + aₙₜ*Kₙₜ*D₁₀ᵀ
    dD₁₀ᵀ = Kₒ*aₒ*D₁₁ᵀ + NT2*aₙₜ*D₀₀ᵀ - O*aₒ*D₁₀ᵀ - aₙₜ*Kₙₜ*D₁₀ᵀ
    dD₀₁ᵀ = -Kₒ*aₒ*D₀₁ᵀ - NT2*aₙₜ*D₀₁ᵀ + O*aₒ*D₀₀ᵀ + aₙₜ*Kₙₜ*D₁₁ᵀ
    dD₁₁ᵀ = -Kₒ*aₒ*D₁₁ᵀ + NT2*aₙₜ*D₀₁ᵀ + O*aₒ*D₁₀ᵀ - aₙₜ*Kₙₜ*D₁₁ᵀ
    dD₀ᴺ = Kₒ*aₒ*D₁ᴺ - O*aₒ*D₀ᴺ
    dD₁ᴺ = -Kₒ*aₒ*D₁ᴺ + O*aₒ*D₀ᴺ
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3



# --- 1. Oct4 overexpress
u00 = zeros(15); u00[2] = 0.5; u00[5] =1; u00[10] =1;u00[14] =1;
p00 = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.]
u0  = SType(u00, p00)
tspan = (0.0,1.5e3)
p = p00

function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p[15] = 0.36# .36
      integrator.p[13] = 0.06 # 0.06
  elseif integrator.t == ts[2]
      integrator.p[15] = 0.0
      integrator.p[13] = 0.0
  end
end
cb = DiscreteCallback(condition, affect!, save_positions=(true,true));

# Solving
prob = ODEProblem(paper_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),callback=cb, tstops=ts)
oct4_over = plot(sol, vars =[:N,:T,:O], color = [:orange :green :blue], title = "CRN Model", ylabel= "Concentration", w = 1.5)
# oct4_over = plot(sol, vars = [5,6,7,8],  title = "CRN Model", ylabel= "Concentration", w = 1.5)
savefig(oct4_over,"/Users/chentianchi/Desktop/epi_plots/oct4_over.png")


# --- 2. TET overexpress
u00 = zeros(15); u00[2] = 0.5; u00[5] =1; u00[10] =1;u00[14] =1;
p00 = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.]
u0  = SType(u00, p00)
tspan = (0.0,1.5e3)
p = p00

function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p[14] = .36
      integrator.p[13] = 0.06
      integrator.p[15] = 0.06
  elseif integrator.t == ts[2]
      integrator.p[14] = 0.0
      integrator.p[13] = 0.0
      integrator.p[15] = 0.0
  end
end
cb = DiscreteCallback(condition, affect!, save_positions=(true,true));

# Solving
prob = ODEProblem(paper_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),callback=cb, tstops=ts)
tet_over = plot(sol, vars = [1,2,9],  title = "CRN Model", ylabel= "Concentration", w = 1.5)
savefig(tet_over,"/Users/chentianchi/Desktop/epi_plots/tet_over.png")





# --- 3. Nanog overexpress
u00 = zeros(15); u00[2] = 0.5; u00[5] =1; u00[10] =1;u00[14] =1;
p00 = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.]
u0  = SType(u00, p00)
tspan = (0.0,1.5e3)
p = p00

function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p[13] = 0.36
      # integrator.p[14] = 0.06
      integrator.p[15] = 0.06
  elseif integrator.t == ts[2]
      integrator.p[13] = 0.0
      # integrator.p[14] = 0.0
      integrator.p[15] = 0.0
  end
end
cb = DiscreteCallback(condition, affect!, save_positions=(true,true));

# Solving
prob = ODEProblem(paper_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),callback=cb, tstops=ts)
nanog_over = plot(sol, vars =[:N,:T,:O], color = [:orange :green :blue], title = "CRN Model", ylabel= "Concentration", w = 1.5)
ylims!(0,1)
savefig(nanog_over,"/Users/chentianchi/Desktop/epi_plots/nanog_over.png")




## steady states
sssol = solve(prob,DynamicSS(Rodas5()), callback=cb, tstops=ts, maxiters = 1e6)
sssol = solve(prob,SSRootfind(), callback=cb, tstops=ts, maxiters = 1e6)
plot(sssol,vars = [1,2,9], legend = true, line=(1.2,:dash))
































# ==================  Paper ODE model
ts = [5e2,7e2]
hill(x,k,n) = x^n/(x^n + k^n)
p_equation = @ode_def begin
    dN = lif + hill(O,KO,1) -N + m1
    dO = lif + hill(O,KO,1)*hill((Kd + N + T)/2 - sqrt(((Kd + N + T)/2)^2 - N*T),Knt,2) -O + m2
    dT = hill(O,KO,1)*hill((Kd + N + T)/2 - sqrt(((Kd + N + T)/2)^2 - N*T),Knt,2) -T + m3
end KO Knt Kd m1 m2 m3 lif

tspan = (0.0,1.5e3)




# ---- oct4 overexpress -----
u0 = SType([0., 0., 0.05], [0.3,0.2,0.1, 0., 0., 0.05, 0.])
p = [0.3,0.2,0.1, 0., 0., 0.05, 0.]

condition2(u,t,integrator) = t in ts
function affect2!(integrator)
  if integrator.t == ts[1]
      integrator.p[5] = .3
      integrator.p[7] = 0.06
  elseif integrator.t == ts[2]
      integrator.p[5] = 0.0
      integrator.p[7] = 0.0
  end
end
cb2 = DiscreteCallback(condition2, affect2!, save_positions=(true,true));

prob = ODEProblem(p_equation,u0,tspan,p)
sol = solve(prob, Tsit5(),callback=cb2, tstops=ts )
# sol = solve(prob,DynamicSS(Rodas5()), callback=cb2, tstops=ts )
oct4_over_p = plot(sol, title = "Paper Model", ylabel= "Concentration", w =1.5)
savefig(oct4_over_p, "/Users/chentianchi/Desktop/epi_plots/oct4_over_p.png")




# ---- nanog overexpress -----
u0 = SType([0., 0., 0.05], [0.3,0.2,0.1, 0., 0., 0.05, 0.])
p = [0.3,0.2,0.1, 0., 0., 0.05, 0.]

condition2(u,t,integrator) = t in ts
function affect2!(integrator)
  if integrator.t == ts[1]
      integrator.p[4] = .3
      integrator.p[7] = 0.06
  elseif integrator.t == ts[2]
      integrator.p[4] = 0.0
      integrator.p[7] = 0.0
  end
end
cb2 = DiscreteCallback(condition2, affect2!, save_positions=(true,true));

prob = ODEProblem(p_equation,u0,tspan,p)
sol = solve(prob, Tsit5(),callback=cb2, tstops=ts )
# sol = solve(prob,DynamicSS(Rodas5()), callback=cb2, tstops=ts )
nanog_over_p = plot(sol, vars =[:N,:T,:O], color = [:orange :green :blue], title = "Paper Model", ylabel= "Concentration", w = 1.5)
ylims!(0,1)
savefig(nanog_over_p, "/Users/chentianchi/Desktop/epi_plots/nanog_over_p.png")


# ---- TET overexpress -----
u0 = SType([0., 0., 0.05], [0.3,0.2,0.1, 0., 0., 0.05, 0.])
p = [0.3,0.2,0.1, 0., 0., 0.05, 0.]

condition2(u,t,integrator) = t in ts
function affect2!(integrator)
  if integrator.t == ts[1]
      integrator.p[6] = .3
      integrator.p[7] = 0.06
  elseif integrator.t == ts[2]
      integrator.p[6] = 0.05
      integrator.p[7] = 0.0
  end
end
cb2 = DiscreteCallback(condition2, affect2!, save_positions=(true,true));

prob = ODEProblem(p_equation,u0,tspan,p)
sol = solve(prob, Tsit5(),callback=cb2, tstops=ts )
# sol = solve(prob,DynamicSS(Rodas5()), callback=cb2, tstops=ts )
TET_over_p = plot(sol, title = "Paper Model", ylabel= "Concentration", w = 1.5)
savefig(TET_over_p, "/Users/chentianchi/Desktop/epi_plots/TET_over_p.png")


# --- Save comparison plots (CRN vs Paper)
TET_compare = plot(tet_over,TET_over_p, layout =2)
savefig(TET_compare, "/Users/chentianchi/Desktop/epi_plots/TET_compare.png")

Oct4_compare = plot(oct4_over,oct4_over_p, layout =2)
savefig(Oct4_compare, "/Users/chentianchi/Desktop/epi_plots/Oct4_compare.png")

Nanog_compare = plot(nanog_over,nanog_over_p, layout =2)
savefig(Nanog_compare, "/Users/chentianchi/Desktop/epi_plots/Nanog_compare.png")
