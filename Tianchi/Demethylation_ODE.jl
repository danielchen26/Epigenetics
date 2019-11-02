using DifferentialEquations, ParameterizedFunctions
using DataFrames, LinearAlgebra, Plots
using Queryverse




# ======= For ODE demethylation Full model =============

# Test Demethy_ode
Demethy_ode = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN  = Dᴺ*α₁ - N*δₒ - N*T*a₁ - N*r₁*Dₙ + NT*Kd*a₁ + r₁*Kₒ*Dᴺ + m1
    dT  = Dᵗ*α₁ - T*δₒ - N*T*a₁ + NT*Kd*a₁ + m2
    dNT = N*T*a₁ - NT*Kd*a₁
    dO  = Dᵒ*α₁ - O*δₒ + Kₒ*Dᴺ*ζ₁ - O*Dₙ*ζ₁ - O*r₁*Dₒ - O*r₁*Dₜ + r₁*Kₒ*Dᵒ + r₁*Kₒ*Dᵗ + m3
    dDₒ = -a0*Dₒ + β*Dₕ - O*r₁*Dₒ + r₁*Kₒ*Dᵒ
    dDᵒ = O*r₁*Dₒ - r₁*Kₒ*Dᵒ
    dD̄ = a0*Dₒ - NT*D̄*k₃
    dDₕ = -β*Dₕ + NT*D̄*k₃
    dDₜ = -NT*k₄*Dₜ - O*r₁*Dₜ + r₁*Kₒ*Dᵗ
    dDᵗ = NT*k₄*Dₜ + O*r₁*Dₜ - r₁*Kₒ*Dᵗ
    dDₙ = Kₒ*Dᴺ*ζ₁ - N*r₁*Dₙ - O*Dₙ*ζ₁ + r₁*Kₒ*Dᴺ
    dDᴺ = -Kₒ*Dᴺ*ζ₁ + N*r₁*Dₙ + O*Dₙ*ζ₁ - r₁*Kₒ*Dᴺ
end Kₒ Kd r₁ a₁ α₁ δₒ a0 β ζ₁ k₃ k₄ m1 m2 m3

p = [0.3, 0.1, 1, 1, 1, 1, 0.5, 1, 1, 1, 1, 0.,0.05,0.]
u0 = rand(12); tspan = (0.,1000.)
prob = ODEProblem(Demethy_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),abstol=1e-8)

plot(sol, vars=[:N, :T,:O], ylims = (0,1.5))



# ======= For ODE demethylation reduced model =============
Demethy_reduced_ode = @ode_def_bare begin
    dN = (N*α₁*r₁ + O*α₁*ζ₁ + Kₒ*m1*r₁ + N*m1*r₁ + Kₒ*m1*ζ₁ + O*m1*ζ₁ - N^2*δₒ*r₁ - Kₒ*N*δₒ*r₁ - Kₒ*N*δₒ*ζ₁ - N*O*δₒ*ζ₁ - Kₒ*N*r₁*ζ₁ + Kₒ*O*r₁*ζ₁)/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)
    dT = m2 - T*δₒ + (α₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)
    dO = (Kd*(Kₒ*a0*β*m3 - Kₒ*O*a0*β*δₒ - Kₒ*O*a0*β*r₁ - Kₒ*O*a0*β*ζ₁ + (Kₒ^2*a0*β*ζ₁*(N*r₁ + O*ζ₁))/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁) + (Kₒ^2*a0*β*r₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄) + (Kₒ*O*a0*β*ζ₁*(N*r₁ + O*ζ₁))/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁) + (N*O*T*α₁*β*k₃)/Kd + (Kₒ*N*T*a0*k₃*m3)/Kd + (Kₒ*N*T*β*k₃*m3)/Kd + (N*O*T*β*k₃*m3)/Kd + (Kₒ*O*a0*β*r₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄) - (N*O^2*T*β*δₒ*k₃)/Kd - (N*O^2*T*β*k₃*r₁)/Kd - (N*O^2*T*β*k₃*ζ₁)/Kd - (Kₒ*N*O*T*a0*δₒ*k₃)/Kd - (Kₒ*N*O*T*β*δₒ*k₃)/Kd - (Kₒ*N*O*T*a0*k₃*r₁)/Kd - (Kₒ*N*O*T*β*k₃*r₁)/Kd - (Kₒ*N*O*T*a0*k₃*ζ₁)/Kd - (Kₒ*N*O*T*β*k₃*ζ₁)/Kd + (Kₒ^2*N*T*a0*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ^2*N*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (N*O^2*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ^2*N*T*a0*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (Kₒ^2*N*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (N*O^2*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (Kₒ*N*O*T*a0*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (2*Kₒ*N*O*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ*N*O*T*a0*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (2*Kₒ*N*O*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄))))/(Kₒ*Kd*a0*β + Kₒ*N*T*a0*k₃ + Kₒ*N*T*β*k₃ + N*O*T*β*k₃)
end Kₒ Kd r₁ α₁ δₒ a0 β ζ₁ k₃ k₄ m1 m2 m3

u0 = rand(3); tspan = (0.,1000.)
p = [0.3, 0.1, 1., 1., 1., 0.5, 1., 1., 1., 1., 0.,0.,0.]
prob = ODEProblem(Demethy_reduced_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),abstol=1e-8)
plot(sol, vars=[:N, :T,:O], ylims = (0,1.5))


anim = @animate for i = 1:1000
    u0 = rand(3)
    prob = ODEProblem(Demethy_reduced_ode,u0,tspan,p)
    sol = solve(prob, Rosenbrock23(),abstol=1e-8)
    plot(sol, vars=[:N, :T,:O], ylims = (0,1.5))
end

gif(anim, "/Users/chentianchi/Desktop/demethy_rd.gif", fps = 10)















# =========================== Below syntax is Maltab Compatible ================================
# =========================== Below syntax is Maltab Compatible ================================
# =========================== Below syntax is Maltab Compatible ================================
# =========================== Below syntax is Maltab Compatible ================================









using DifferentialEquations,  ParameterizedFunctions
using Plots

module epi
    include("functions.jl")
    include("./Model_collection/ODE_model_set.jl")
end


# ======= For ODE demethylation Full model =============
Demethy_ode_MatlabC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN    = DN*alpha1 - N*delta0 - N*T*a1 - N*r1*Dn + NT*Kd*a1 + r1*K0*DN + m1
    dT    = DT*alpha1 - T*delta0 - N*T*a1 + NT*Kd*a1 + m2
    dNT   = N*T*a1 - NT*Kd*a1
    dO    = DO*alpha1 - O*delta0 + K0*DN*xi1 - O*Dn*xi1 - O*r1*Do - O*r1*Dt + r1*K0*DO + r1*K0*DT + m3
    dDo   = -a0*Do + beta*D5hmc - O*r1*Do + r1*K0*DO
    dDO   = O*r1*Do - r1*K0*DO
    dD5mc = a0*Do - k3*NT*D5mc
    dD5hmc= -beta*D5hmc + k3*NT*D5mc
    dDt   = -O*r1*Dt + r1*K0*DT
    dDT   = O*r1*Dt - r1*K0*DT
    dDn   = K0*DN*xi1 - N*r1*Dn - O*Dn*xi1 + r1*K0*DN
    dDN   = -K0*DN*xi1 + N*r1*Dn + O*Dn*xi1 - r1*K0*DN
end K0 Kd r1 a1 alpha1 delta0 a0 beta xi1 k3 m1 m2 m3




p = [0.3, 0.1, 1, 1, 1, 1, 0.5, 1, 1, 1, 0.,0.05,0.]
tspan = (0.,1000.)

N_sample = 10
anim = @animate for u01 in epi.randfixsum(N_sample,4,1) , u02 in epi.randfixsum(N_sample,2,1), u03 in epi.randfixsum(N_sample,2,1)
    u0 = [rand(1,4) u01 u02 u03]
    prob = ODEProblem(Demethy_ode_MatlabC, u0, tspan, p)
    sol  = solve(prob, Rosenbrock23())
    plot(sol, vars =[:N, :T, :O], ylims = (0,2))
end
gif(anim, fps = 10)




# ======= For ODE demethylation reduced model =============
Demethy_reduced_ode_MatlabC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN = (N*alpha1*r1 + O*alpha1*xi1 + K0*m1*r1 + N*m1*r1 + K0*m1*xi1 + O*m1*xi1 - N^2*delta0*r1 - K0*N*delta0*r1 - K0*N*delta0*xi1 - N*O*delta0*xi1 - K0*N*r1*xi1 + K0*O*r1*xi1)/(K0*r1 + N*r1 + K0*xi1 + O*xi1)

    dT = m2 - T*delta0 + (O*alpha1)/(K0 + O)

    dO = (Kd*(K0*a0*beta*m3 - K0*O*a0*beta*delta0 - K0*O*a0*beta*xi1 + (K0^2*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (K0*O*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (N*O*T*alpha1*beta*k3)/Kd + (K0*N*T*a0*k3*m3)/Kd + (K0*N*T*beta*k3*m3)/Kd + (N*O*T*beta*k3*m3)/Kd - (N*O^2*T*beta*delta0*k3)/Kd - (N*O^2*T*beta*k3*xi1)/Kd - (K0*N*O*T*a0*delta0*k3)/Kd - (K0*N*O*T*beta*delta0*k3)/Kd - (K0*N*O*T*a0*k3*xi1)/Kd - (K0*N*O*T*beta*k3*xi1)/Kd + (K0^2*N*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0^2*N*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (N*O^2*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0*N*O*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (2*K0*N*O*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1))))/(K0*Kd*a0*beta + K0*N*T*a0*k3 + K0*N*T*beta*k3 + N*O*T*beta*k3)
end K0 Kd r1 alpha1 delta0 a0 beta xi1 k3 m1 m2 m3


u0 = rand(3); tspan = (0.,1000.)
# p = [0.3, 0.1, 1., 1., 0.5, 1., 1., 1., 1., 0.,0.05,0.]
p = [0.3, 0.1, 1.09178, 2.86155, 2.86159, 1.28656, 2.85916, 2.53682, 2.56739, 0.,0.05,0.]
prob = ODEProblem(Demethy_reduced_ode_MatlabC,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),abstol=1e-8)
plot(sol, vars=[:N, :T,:O], ylims = (0,2.5))
sol[end]

anim = @animate for i = 1:1000
    u0 = rand(3)
    prob = ODEProblem(Demethy_reduced_ode_MatlabC,u0,tspan,p)
    sol = solve(prob, Rosenbrock23(),abstol=1e-8)
    plot(sol, vars=[:N, :T,:O], ylims = (0,2.5))
end

gif(anim, "/Users/chentianchi/Desktop/demethy_rd_NO_NT_TET.gif", fps = 5)
