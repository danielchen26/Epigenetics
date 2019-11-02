# ODE model set
# --- Collection of all ODE models converted from CRN
using DifferentialEquations, ParameterizedFunctions


# 1. Paper's model
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


# 2. Paper's reduced model
paper_reduce_ode = @ode_def begin

    dN = -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1

    dT = (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))

    dO = (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3




# 3. Methylation model | with last rn : NT activate TET
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





## 4. Methylation reduced model  | with last rn : NT activate TET
Demethy_reduced_ode = @ode_def_bare begin
    dN = (N*α₁*r₁ + O*α₁*ζ₁ + Kₒ*m1*r₁ + N*m1*r₁ + Kₒ*m1*ζ₁ + O*m1*ζ₁ - N^2*δₒ*r₁ - Kₒ*N*δₒ*r₁ - Kₒ*N*δₒ*ζ₁ - N*O*δₒ*ζ₁ - Kₒ*N*r₁*ζ₁ + Kₒ*O*r₁*ζ₁)/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)

    dT = m2 - T*δₒ + (α₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)

    dO = (Kd*(Kₒ*a0*β*m3 - Kₒ*O*a0*β*δₒ - Kₒ*O*a0*β*r₁ - Kₒ*O*a0*β*ζ₁ + (Kₒ^2*a0*β*ζ₁*(N*r₁ + O*ζ₁))/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁) + (Kₒ^2*a0*β*r₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄) + (Kₒ*O*a0*β*ζ₁*(N*r₁ + O*ζ₁))/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁) + (N*O*T*α₁*β*k₃)/Kd + (Kₒ*N*T*a0*k₃*m3)/Kd + (Kₒ*N*T*β*k₃*m3)/Kd + (N*O*T*β*k₃*m3)/Kd + (Kₒ*O*a0*β*r₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄) - (N*O^2*T*β*δₒ*k₃)/Kd - (N*O^2*T*β*k₃*r₁)/Kd - (N*O^2*T*β*k₃*ζ₁)/Kd - (Kₒ*N*O*T*a0*δₒ*k₃)/Kd - (Kₒ*N*O*T*β*δₒ*k₃)/Kd - (Kₒ*N*O*T*a0*k₃*r₁)/Kd - (Kₒ*N*O*T*β*k₃*r₁)/Kd - (Kₒ*N*O*T*a0*k₃*ζ₁)/Kd - (Kₒ*N*O*T*β*k₃*ζ₁)/Kd + (Kₒ^2*N*T*a0*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ^2*N*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (N*O^2*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ^2*N*T*a0*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (Kₒ^2*N*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (N*O^2*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (Kₒ*N*O*T*a0*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (2*Kₒ*N*O*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ*N*O*T*a0*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (2*Kₒ*N*O*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄))))/(Kₒ*Kd*a0*β + Kₒ*N*T*a0*k₃ + Kₒ*N*T*β*k₃ + N*O*T*β*k₃)
end Kₒ Kd r₁ α₁ δₒ a0 β ζ₁ k₃ k₄ m1 m2 m3





# =========================== Below syntax is Maltab NC ================================
# =========================== Below syntax is Maltab NC ================================
# =========================== Below syntax is Maltab NC ================================
# =========================== Below syntax is Maltab NC ================================


# 5. Methylation model | without last rn : NT activate TET
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








Demethy_reduced_ode_MatlabC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O

    dN = (N*alpha1*r1 + O*alpha1*xi1 + K0*m1*r1 + N*m1*r1 + K0*m1*xi1 + O*m1*xi1 - N^2*delta0*r1 - K0*N*delta0*r1 - K0*N*delta0*xi1 - N*O*delta0*xi1 - K0*N*r1*xi1 + K0*O*r1*xi1)/(K0*r1 + N*r1 + K0*xi1 + O*xi1)

    dT = m2 - T*delta0 + (O*alpha1)/(K0 + O)

    dO = (Kd*(K0*a0*beta*m3 - K0*O*a0*beta*delta0 - K0*O*a0*beta*xi1 + (K0^2*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (K0*O*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (N*O*T*alpha1*beta*k3)/Kd + (K0*N*T*a0*k3*m3)/Kd + (K0*N*T*beta*k3*m3)/Kd + (N*O*T*beta*k3*m3)/Kd - (N*O^2*T*beta*delta0*k3)/Kd - (N*O^2*T*beta*k3*xi1)/Kd - (K0*N*O*T*a0*delta0*k3)/Kd - (K0*N*O*T*beta*delta0*k3)/Kd - (K0*N*O*T*a0*k3*xi1)/Kd - (K0*N*O*T*beta*k3*xi1)/Kd + (K0^2*N*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0^2*N*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (N*O^2*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0*N*O*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (2*K0*N*O*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1))))/(K0*Kd*a0*beta + K0*N*T*a0*k3 + K0*N*T*beta*k3 + N*O*T*beta*k3)
end K0 Kd r1 alpha1 delta0 a0 beta xi1 k3 m1 m2 m3





# convert from Demethy_TF_MC.f_symfuncs
Demethy_TF_ode_MC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN    = m1 + Dn1*alphaN - N*delta - N*T*a1 + NT*Kd*a1
    dT    = m2 + Dt11*alphaT - T*delta - N*T*a1 + NT*Kd*a1
    dNT   = -d*NT^2 + 2*d*NT2 + N*T*a1 - NT*Kd*a1
    dNT2  = (1/2)*d*NT^2 - d*NT2 + K_nt*a_nt*Do10 + K_nt*a_nt*Do11 + K_nt*a_nt*Dt10 + K_nt*a_nt*Dt11 - NT2*a_nt*Do00 - NT2*a_nt*Do01 - NT2*a_nt*Dt00 - NT2*a_nt*Dt01
    dDo00 = -a_dn*Do00 + beta*D5hmc + KO*aO*Do01 + K_nt*a_nt*Do10 - NT2*a_nt*Do00 - O*aO*Do00
    dDo10 = KO*aO*Do11 - K_nt*a_nt*Do10 + NT2*a_nt*Do00 - O*aO*Do10
    dDo01 = -KO*aO*Do01 + K_nt*a_nt*Do11 - NT2*a_nt*Do01 + O*aO*Do00
    dDo11 = -KO*aO*Do11 - K_nt*a_nt*Do11 + NT2*a_nt*Do01 + O*aO*Do10
    dO    = m3 + Do11*alphaO - O*delta + KO*aO*Dn1 + KO*aO*Do01 + KO*aO*Do11 + KO*aO*Dt01 + KO*aO*Dt11 - O*aO*Dn0 - O*aO*Do00 - O*aO*Do10 - O*aO*Dt00 - O*aO*Dt10
    dDt00 = KO*aO*Dt01 + K_nt*a_nt*Dt10 - NT2*a_nt*Dt00 - O*aO*Dt00
    dDt10 = KO*aO*Dt11 - K_nt*a_nt*Dt10 + NT2*a_nt*Dt00 - O*aO*Dt10
    dDt01 = -KO*aO*Dt01 + K_nt*a_nt*Dt11 - NT2*a_nt*Dt01 + O*aO*Dt00
    dDt11 = -KO*aO*Dt11 - K_nt*a_nt*Dt11 + NT2*a_nt*Dt01 + O*aO*Dt10
    dDn0  = KO*aO*Dn1 - O*aO*Dn0
    dDn1  = -KO*aO*Dn1 + O*aO*Dn0
    dD5mc = a_dn*Do00 - kh*NT*D5mc
    dD5hmc = -beta*D5hmc + kh*NT*D5mc
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn beta kh m1 m2 m3




# Demethy + TF reduced to 3d(NTO)
Demethy_TF_ode_reduced_MC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN = m1 - N*delta + (O*alphaN)/(KO + O)
    dT = m2 - T*delta + (N^2*O*T^2*alphaT)/((2*K_nt*Kd^2 + N^2*T^2)*(KO + O))
    dO = (N^3*O*T^3*alphaO*beta*kh + KO*N^3*T^3*beta*kh*m3 + N^3*O*T^3*beta*kh*m3 - N^3*O^2*T^3*beta*delta*kh + 2*KO*K_nt*Kd^3*a_dn*beta*m3 - 2*KO*K_nt*Kd^3*O*a_dn*beta*delta - KO*N^3*O*T^3*beta*delta*kh + 2*KO*K_nt*Kd^2*N*T*a_dn*kh*m3 + 2*KO*K_nt*Kd^2*N*T*beta*kh*m3 + 2*K_nt*Kd^2*N*O*T*beta*kh*m3 - 2*K_nt*Kd^2*N*O^2*T*beta*delta*kh - 2*KO*K_nt*Kd^2*N*O*T*a_dn*delta*kh - 2*KO*K_nt*Kd^2*N*O*T*beta*delta*kh)/(2*KO*K_nt*Kd^3*a_dn*beta + KO*N^3*T^3*beta*kh + N^3*O*T^3*beta*kh + 2*KO*K_nt*Kd^2*N*T*a_dn*kh + 2*KO*K_nt*Kd^2*N*T*beta*kh + 2*K_nt*Kd^2*N*O*T*beta*kh)
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn beta kh m1 m2 m3
