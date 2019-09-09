using DifferentialEquations
using ParameterizedFunctions
using SymEngine

N, T, NT, NT2, D₀₀ᴼ, D₁₀ᴼ, D₀₁ᴼ, D₁₁ᴼ, O, D₀₀ᵀ, D₁₀ᵀ, D₀₁ᵀ, D₁₁ᵀ, D₀ᴺ, D₁ᴺ = @vars N T NT NT2 D₀₀ᴼ D₁₀ᴼ D₀₁ᴼ D₁₁ᴼ O D₀₀ᵀ D₁₀ᵀ D₀₁ᵀ D₁₁ᵀ D₀ᴺ D₁ᴺ

u0 = [N; T; NT; NT2; D₀₀ᴼ; D₁₀ᴼ; D₀₁ᴼ; D₁₁ᴼ; O; D₀₀ᵀ; D₁₀ᵀ; D₀₁ᵀ; D₁₁ᵀ; D₀ᴺ; D₁ᴺ]

crn_ode = @ode_def begin
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

p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.]
prob = ODEProblem(crn_ode,u0,(0.0,1.0), p)
sol = solve(prob,RK4(),dt=1/2,adaptive=false)
end_solution  = lambdify(sol[end])


#  ------ Simplify equation  ----- paper model ------
using SymPy

@vars N T NT NT2 D₀₀ᴼ D₁₀ᴼ D₀₁ᴼ D₁₁ᴼ O D₀₀ᵀ D₁₀ᵀ D₀₁ᵀ D₁₁ᵀ D₀ᴺ D₁ᴺ nonnegative=true a real=true
Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3= 0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.


fN   = -N*δ + αₙ*D₁ᴺ - N*T*a₁ + NT*Kd*a₁ + m1
fT   = -T*δ + αₜ*D₁₁ᵀ - N*T*a₁ + NT*Kd*a₁ + m2
fO   = -O*δ + αₒ*D₁₁ᴼ + Kₒ*aₒ*D₀₁ᴼ + Kₒ*aₒ*D₀₁ᵀ + Kₒ*aₒ*D₁ᴺ + Kₒ*aₒ*D₁₁ᴼ + Kₒ*aₒ*D₁₁ᵀ - O*aₒ*D₀ᴺ - O*aₒ*D₀₀ᴼ - O*aₒ*D₀₀ᵀ - O*aₒ*D₁₀ᴼ - O*aₒ*D₁₀ᵀ + m3
fD₀₀ᴼ = Kₒ*aₒ*D₀₁ᴼ - NT2*aₙₜ*D₀₀ᴼ - O*aₒ*D₀₀ᴼ + aₙₜ*Kₙₜ*D₁₀ᴼ
fD₁₀ᴼ = Kₒ*aₒ*D₁₁ᴼ + NT2*aₙₜ*D₀₀ᴼ - O*aₒ*D₁₀ᴼ - aₙₜ*Kₙₜ*D₁₀ᴼ
fD₀₁ᴼ = -Kₒ*aₒ*D₀₁ᴼ - NT2*aₙₜ*D₀₁ᴼ + O*aₒ*D₀₀ᴼ + aₙₜ*Kₙₜ*D₁₁ᴼ
fD₁₁ᴼ = -Kₒ*aₒ*D₁₁ᴼ + NT2*aₙₜ*D₀₁ᴼ + O*aₒ*D₁₀ᴼ - aₙₜ*Kₙₜ*D₁₁ᴼ
fD₀₀ᵀ = Kₒ*aₒ*D₀₁ᵀ - NT2*aₙₜ*D₀₀ᵀ - O*aₒ*D₀₀ᵀ + aₙₜ*Kₙₜ*D₁₀ᵀ
fD₁₀ᵀ = Kₒ*aₒ*D₁₁ᵀ + NT2*aₙₜ*D₀₀ᵀ - O*aₒ*D₁₀ᵀ - aₙₜ*Kₙₜ*D₁₀ᵀ
fD₀₁ᵀ = -Kₒ*aₒ*D₀₁ᵀ - NT2*aₙₜ*D₀₁ᵀ + O*aₒ*D₀₀ᵀ + aₙₜ*Kₙₜ*D₁₁ᵀ
fD₁₁ᵀ = -Kₒ*aₒ*D₁₁ᵀ + NT2*aₙₜ*D₀₁ᵀ + O*aₒ*D₁₀ᵀ - aₙₜ*Kₙₜ*D₁₁ᵀ
fD₀ᴺ  = Kₒ*aₒ*D₁ᴺ - O*aₒ*D₀ᴺ
fD₁ᴺ  = -Kₒ*aₒ*D₁ᴺ + O*aₒ*D₀ᴺ
fNT  = -NT*β - d*NT^2 + 2*d*NT2 + N*T*a₁ - NT*Kd*a₁
fNT2 = (1/2)*d*NT^2 - d*NT2 - β*NT2 - NT2*aₙₜ*D₀₀ᴼ - NT2*aₙₜ*D₀₀ᵀ - NT2*aₙₜ*D₀₁ᴼ - NT2*aₙₜ*D₀₁ᵀ + aₙₜ*Kₙₜ*D₁₀ᴼ + aₙₜ*Kₙₜ*D₁₀ᵀ + aₙₜ*Kₙₜ*D₁₁ᴼ + aₙₜ*Kₙₜ*D₁₁ᵀ
@show exs = [ fN, fT, fO,  fD₀₀ᴼ, fD₁₀ᴼ, fD₀₁ᴼ, fD₁₁ᴼ, fD₀₀ᵀ, fD₁₀ᵀ, fD₀₁ᵀ, fD₁₁ᵀ, fD₀ᴺ, fD₁ᴺ, fNT, fNT2]





# 1. Substitute with conservation law
D₁₁ᴼ_sub = 1 - D₀₀ᴼ - D₁₀ᴼ - D₀₁ᴼ
D₁₁ᵀ_sub = 1 - D₀₀ᵀ - D₁₀ᵀ - D₀₁ᵀ
D₁ᴺ_sub  = 1 - D₀ᴺ
@show exs = exs.subs.(D₁ᴺ, D₁ᴺ_sub )
@show exs =exs.subs.(D₁₁ᵀ, D₁₁ᵀ_sub )
@show exs =exs.subs.(D₁₁ᴼ, D₁₁ᴼ_sub )
@. exs = simplify(exs)

# 2. Solve for other variables except for 1:=> N, 2:=>T, 3:=>O.
dD₀₀ᴼ = solve(fD₀₀ᴼ, D₀₀ᴼ)
@show exs =exs.subs.(D₀₀ᴼ, dD₀₀ᴼ)[1]
dD₁₀ᴼ = solve(fD₁₀ᴼ, D₁₀ᴼ)
@show exs =exs.subs.(D₁₀ᴼ, dD₁₀ᴼ)[1]
dD₀₁ᴼ = solve(fD₀₁ᴼ, D₀₁ᴼ)
@show exs =exs.subs.(D₀₁ᴼ, dD₀₁ᴼ)[1]
@. exs = simplify(exs)

dD₀₀ᵀ = solve(fD₀₀ᵀ, D₀₀ᵀ)
@show exs =exs.subs.(D₀₀ᵀ, dD₀₀ᵀ)[1]
dD₁₀ᵀ = solve(fD₁₀ᵀ, D₁₀ᵀ)
@show exs =exs.subs.(D₁₀ᵀ, dD₁₀ᵀ)[1]
dD₀₁ᵀ = solve(fD₀₁ᵀ, D₀₁ᵀ)
@show exs =exs.subs.(D₀₁ᵀ, dD₀₁ᵀ)[1]
@. exs = simplify(exs)

dD₀ᴺ = solve(fD₀ᴺ, D₀ᴺ)
@show exs =exs.subs.(D₀ᴺ, dD₀ᴺ)[1]
@. exs = simplify(exs)

dNT = solve(fNT, NT)
@show exs =exs.subs.(NT, dNT)[1]
dNT2 = solve(fNT2, NT2)
@show exs =exs.subs.(NT2, dNT2)[1]
@. exs = simplify(exs)

@show exs







# @show d = solve(exs)
using DataFrames
df = DataFrame(d[1])
@show eachcol(df[:,1:end],true)
