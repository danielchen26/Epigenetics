# ---------------   Testing the ode control example
using DifferentialEquations, ParameterizedFunctions,Plots;gr()
mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::Array{T,1}
end

ts = [5,10]
condition(u,t,integrator) = t in ts
function affect!(integrator)
  if integrator.t == ts[1]
      integrator.p[26] = 1.0
  elseif integrator.t == ts[2]
      integrator.p[26] = 0.0
  end
end

save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions);


O_N_auto_ode = @ode_def_bare begin
    dO = -O*β - O*δₒ - O^2*a₀ + 2*O₂*d₀ + k*Dₒ*α₁ + m1;
    dO₂ = Dᴺ*ζ₀ + (1/2)*O^2*a₀ - O₂*d₀ + rᵦ*D♇ + r₀*Dₒ - β*O₂ - D*r₁*O₂ - O₂*Dₙ*ζ₁ - r₂*O₂*Dₜ;
    dD = -D*aₐ + r₀*Dₒ + β*Dₕ - D*r₁*O₂ + NT*D̄*k₃;
    dDₒ = -r₀*Dₒ + D*r₁*O₂;
    dD̄ = D*aₐ + d₂*C₂ - NT*D̄*k₃ - TET*D̄*a₂;
    dTET = D♇*ϕ₁ + NT*η₀ - TET*δₜ + d₂*C₂ + k₂*C₂ - β*TET - N*TET*η₁ - TET*D̄*a₂ + m3;
    dC₂ = -d₂*C₂ - k₂*C₂ - β*C₂ + TET*D̄*a₂;
    dDₕ = k₂*C₂ - β*Dₕ;
    dDₜ = rᵦ*D♇ - r₂*O₂*Dₜ;
    dD♇ = -rᵦ*D♇ + r₂*O₂*Dₜ;
    dN = -N*β -N*δₒ -N^2*a₀ + NT*η₀ + 2*N₂*d₀ - N*TET*η₁ + k*Dᴺ*α₁ + m2;
    dN₂ = (1/2)*N^2*a₀ - N₂*d₀ + r₀*Dᴺ - β*N₂ - r₁*N₂*Dₙ;
    dDₙ =  Dᴺ*ζ₀ + r₀*Dᴺ - O₂*Dₙ*ζ₁ - r₁*N₂*Dₙ;
    dDᴺ = -Dᴺ*ζ₀ - r₀*Dᴺ + O₂*Dₙ*ζ₁ + r₁*N₂*Dₙ;
    dNT = -NT*η₀ + (N + m2)*TET*η₁
end a₀ d₀ r₁ r₀ α₀ α₁ k  β  δₒ aₐ a₂ d₂ k₂ ϕ₀ ϕ₁ r₂ rᵦ δₜ ζ₁ ζ₀ η₁ η₀ k₃ m1 m2 m3

Array([p2;[0,0,0]])

# Solving
si = zeros(Float64,15)
si[[3,4,6,10,13,14]] .= 0.5
u0 = SType(si, Array([p2;[0,0,0]]))
tspan = (0.0,1e3)
p = Array([p2;[0,0,0]])
prob = ODEProblem(O_N_auto_ode,u0,tspan,p)

sol = solve(prob,Tsit5(),callback=cb,tstops=ts, maxiters = 1e6)
plot(sol, vars = [1,6,11] )

sssol = solve(prob,DynamicSS(Tsit5()))
plot(sssol, vars = getindex(O_N_auto_ode.syms, [1,6,11]),  legend = true, line=(1.2,:dash))




#  Animation
si = zeros(Float64,15)
tspan = (0.0,1e5)
anim = @animate for i=0:0.5:10

    ts = [5,10]
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
      if integrator.t == ts[1]
          integrator.p[26] = i
      elseif integrator.t == ts[2]
          integrator.p[26] = 0.0
      end
    end
    save_positions = (true,true)
    cb = DiscreteCallback(condition, affect!, save_positions=save_positions);

    si[[3,4,6,10,13,14]] .= .5
    u0 = SType(si, Array([p2;[0,0,0]]))
    p = Array([p2;[0,0,0]])
    prob = ODEProblem(O_N_auto_ode,u0,tspan,p)
    # steady states solutions
    sssol = solve(prob,DynamicSS(Tsit5()), callback=cb,tstops=ts)
    plot(sssol, vars = getindex(O_N_auto_ode.syms, [1,6,11]),  legend = true, line=(1.2,:dash))
    # ylims!(0,1)
    title!("Model Evolution : Strength Response : DM_initial = $i ")
    xaxis!("Time")
    yaxis!("Concentration")
end
