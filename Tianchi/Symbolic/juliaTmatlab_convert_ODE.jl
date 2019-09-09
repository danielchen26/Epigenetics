using SymPy

# Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3= 0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.
# all syms
@vars N T NT NT2 D₀₀ᴼ D₁₀ᴼ D₀₁ᴼ D₁₁ᴼ O D₀₀ᵀ D₁₀ᵀ D₀₁ᵀ D₁₁ᵀ D₀ᴺ D₁ᴺ Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3
# ------ paper_reduce_ode
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
@show exs = [ fN, fT, fO, fD₀₀ᴼ, fD₁₀ᴼ, fD₀₁ᴼ, fD₁₁ᴼ, fD₀₀ᵀ, fD₁₀ᵀ, fD₀₁ᵀ, fD₁₁ᵀ, fD₀ᴺ, fD₁ᴺ, fNT, fNT2]


# Convert variables
@vars Do00 Do10 Do01 Do11 Dt00 Dt10 Dt01 Dt11 Dn0 Dn1
@show exs = exs.subs(D₀₀ᴼ, Do00)
@show exs = exs.subs(D₁₀ᴼ, Do10)
@show exs = exs.subs(D₀₁ᴼ, Do01)
@show exs = exs.subs(D₁₁ᴼ, Do11)

@show exs = exs.subs(D₀₀ᵀ, Dt00)
@show exs = exs.subs(D₁₀ᵀ, Dt10)
@show exs = exs.subs(D₀₁ᵀ, Dt01)
@show exs = exs.subs(D₁₁ᵀ, Dt11)

@show exs = exs.subs(D₀ᴺ, Dn0)
@show exs = exs.subs(D₁ᴺ, Dn1)

# Convert parameters
@vars K_o K_nt Kd a1 d a_nt a_0 alpha_t alpha_o alpha_n delta beta
@show exs = exs.subs(Kₒ, K_o)
@show exs = exs.subs(Kₙₜ, K_nt)
@show exs = exs.subs(Kd, Kd)
@show exs = exs.subs(a₁, a1)
@show exs = exs.subs(d, d)
@show exs = exs.subs(aₙₜ, a_nt)
@show exs = exs.subs(aₒ, a_0)
@show exs = exs.subs(αₜ, alpha_t)
@show exs = exs.subs(αₒ, alpha_o)
@show exs = exs.subs(αₙ, alpha_n)
@show exs = exs.subs(δ, delta)
@show exs = exs.subs(β, beta)


display(exs)





# paper_reduce_ode
@vars N T O  K_o  K_nt  Kd  a1  d  a_nt  a_0  alpha_t  alpha_o alpha_n delta beta m1 m2 m3
dN = m1 - N*delta - N*T*a1 + (O*alpha_n)/(K_o + O) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d)

dT = -(K_nt*K_o*T*delta - K_nt*O*m2 - K_nt*K_o*m2 + K_nt*O*T*delta + K_nt*K_o*N*T*a1 + K_nt*N*O*T*a1 - (O*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (O*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*O*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*T*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (O*T*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*O*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*N*T*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (N*O*T*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + O))

dO = (K_nt*K_o*m3 + K_nt*O*m3 - K_nt*O^2*delta - K_nt*K_o*O*delta + (O*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (O*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (O^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*O*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + O))


DR = [dN, dT, dO]



# Convert parameters
@vars N T O  Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3
@show DR = DR.subs(K_o , Kₒ )
@show DR = DR.subs(K_nt, Kₙₜ )
@show DR = DR.subs(Kd , Kd)
@show DR = DR.subs(a1,a₁)
@show DR = DR.subs(d ,d )
@show DR = DR.subs(a_nt ,aₙₜ)
@show DR = DR.subs(a_0  ,aₒ)
@show DR = DR.subs(alpha_t,αₜ)
@show DR = DR.subs(alpha_o,αₒ)
@show DR = DR.subs(alpha_n,αₙ)
@show DR = DR.subs(delta,δ)
@show DR = DR.subs(beta ,β)

display(DR)




# for matlab
dx(1) = m1 - x(1)*delta - x(1)*x(2)*a1 + (x(3)*alpha_n)/(K_o + x(3)) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d)

dx(2) = -(K_nt*K_o*x(2)*delta - K_nt*x(3)*m2 - K_nt*K_o*m2 + K_nt*x(3)*x(2)*delta + K_nt*K_o*x(1)*x(2)*a1 + K_nt*x(1)*x(3)*x(2)*a1 - (x(3)*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*x(3)*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*x(3)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*x(1)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(1)*x(3)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)))

dx(3) = (K_nt*K_o*m3 + K_nt*x(3)*m3 - K_nt*x(3)^2*delta - K_nt*K_o*x(3)*delta + (x(3)*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*x(3)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)))








# numerator of f1f2f3
@vars N T O  Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3  K_o  K_nt  Kd  a1  d  a_nt  a_0  alpha_t  alpha_o alpha_n delta beta
fN = 2*O*alpha_n*beta*d - Kd*O*a1*beta^2 - K_o*Kd^2*a1^2*beta - K_o*Kd^2*a1^2*d - Kd^2*O*a1^2*beta - Kd^2*O*a1^2*d - K_o*Kd*a1*beta^2 + 2*K_o*beta*d*m1 + 2*O*beta*d*m1 - K_o*Kd*a1*beta*d - Kd*O*a1*beta*d - 2*K_o*N*beta*d*delta - 2*N*O*beta*d*delta + K_o*Kd*a1*beta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + K_o*Kd*a1*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + Kd*O*a1*beta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + Kd*O*a1*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 2*K_o*N*T*a1*beta*d - 2*N*O*T*a1*beta*d
fT = - 4*K_o*Kd*a1*beta^5 - 4*Kd*O*a1*beta^5 + 4*O*alpha_t*beta^4*d + 4*K_o*beta^4*d*m2 + 4*O*beta^4*d*m2 + 4*O*alpha_t*beta^3*d^2 + 4*K_o*beta^3*d^2*m2 + 4*O*beta^3*d^2*m2 - 12*K_o*Kd^2*a1^2*beta^4 - 12*K_o*Kd^3*a1^3*beta^3 - 4*K_o*Kd^4*a1^4*beta^2 - 4*K_o*Kd^4*a1^4*d^2 - 12*Kd^2*O*a1^2*beta^4 - 12*Kd^3*O*a1^3*beta^3 - 4*Kd^4*O*a1^4*beta^2 - 4*Kd^4*O*a1^4*d^2 + 6*K_o*Kd^2*a1^2*beta^3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 3*K_o*Kd^3*a1^3*beta^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 3*K_o*Kd^3*a1^3*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 6*Kd^2*O*a1^2*beta^3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 3*Kd^3*O*a1^3*beta^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 3*Kd^3*O*a1^3*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 24*K_o*Kd^2*a1^2*beta^3*d - 12*K_o*Kd^3*a1^3*beta*d^2 - 24*K_o*Kd^3*a1^3*beta^2*d - 24*Kd^2*O*a1^2*beta^3*d - 12*Kd^3*O*a1^3*beta*d^2 - 24*Kd^3*O*a1^3*beta^2*d + K_o*Kd*a1*beta^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(3/2) + 3*K_o*Kd*a1*beta^4*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + K_o*Kd*a1*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(3/2) + Kd*O*a1*beta^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(3/2) + 3*Kd*O*a1*beta^4*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + Kd*O*a1*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(3/2) - 4*O*alpha_t*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*K_o*beta^3*d*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*O*beta^3*d*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 8*K_o*Kd*a1*beta^4*d - 8*Kd*O*a1*beta^4*d - 4*K_o*T*beta^4*d*delta - 4*O*T*beta^4*d*delta - 12*K_o*Kd^2*a1^2*beta^2*d^2 - 12*Kd^2*O*a1^2*beta^2*d^2 - 4*O*alpha_t*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*K_o*beta^2*d^2*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*O*beta^2*d^2*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*K_o*Kd*a1*beta^3*d^2 - 8*K_o*Kd^4*a1^4*beta*d - 4*Kd*O*a1*beta^3*d^2 - 8*Kd^4*O*a1^4*beta*d + 16*K_nt*K_o*beta^3*d^2*m2 - 4*K_o*T*beta^3*d^2*delta + 16*K_nt*O*beta^3*d^2*m2 - 4*O*T*beta^3*d^2*delta + 3*K_o*Kd*a1*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 6*K_o*Kd^3*a1^3*beta*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 3*Kd*O*a1*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 6*Kd^3*O*a1^3*beta*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*T*beta^2*d^2*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*O*T*beta^2*d^2*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 8*K_nt*K_o*Kd*a1*beta^3*d^2 - 8*K_nt*Kd*O*a1*beta^3*d^2 - 4*K_o*N*T*a1*beta^3*d^2 - 16*K_nt*K_o*T*beta^3*d^2*delta - 4*N*O*T*a1*beta^3*d^2 - 16*K_nt*O*T*beta^3*d^2*delta + 8*Kd*O*a1*alpha_t*beta^2*d^2 + 8*K_o*Kd*a1*beta^2*d^2*m2 + 2*K_o*Kd*a1*beta*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(3/2) + 8*Kd*O*a1*beta^2*d^2*m2 + 2*Kd*O*a1*beta*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(3/2) - 8*K_o*N^2*T^2*a1^2*beta^2*d^2 - 8*N^2*O*T^2*a1^2*beta^2*d^2 + 6*K_o*Kd^2*a1^2*beta*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 12*K_o*Kd^2*a1^2*beta^2*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 6*Kd^2*O*a1^2*beta*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 12*Kd^2*O*a1^2*beta^2*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 8*K_nt*K_o*Kd^2*a1^2*beta^3*d - 8*K_nt*Kd^2*O*a1^2*beta^3*d + 4*Kd^2*O*a1^2*alpha_t*beta*d^2 + 4*Kd^2*O*a1^2*alpha_t*beta^2*d + 4*K_o*Kd^2*a1^2*beta*d^2*m2 + 4*K_o*Kd^2*a1^2*beta^2*d*m2 + 6*K_o*Kd*a1*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*Kd^2*O*a1^2*beta*d^2*m2 + 4*Kd^2*O*a1^2*beta^2*d*m2 + 6*Kd*O*a1*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*T*beta^3*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*O*T*beta^3*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 8*K_nt*K_o*Kd*a1*beta^4*d - 8*K_nt*Kd*O*a1*beta^4*d - 4*K_o*N*T*a1*beta^4*d - 4*N*O*T*a1*beta^4*d + 8*Kd*O*a1*alpha_t*beta^3*d - 8*K_nt*K_o*Kd^2*a1^2*beta^2*d^2 + 8*K_o*Kd*a1*beta^3*d*m2 - 8*K_nt*Kd^2*O*a1^2*beta^2*d^2 + 8*Kd*O*a1*beta^3*d*m2 + 4*K_o*N*T*a1*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*N*O*T*a1*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*Kd*O*a1*alpha_t*beta*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*Kd*O*a1*alpha_t*beta^2*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*K_o*Kd*a1*beta*d^2*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*K_o*Kd*a1*beta^2*d*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*Kd*O*a1*beta*d^2*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*Kd*O*a1*beta^2*d*m2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 8*K_o*Kd*T*a1*beta^3*d*delta - 8*Kd*O*T*a1*beta^3*d*delta + 8*K_nt*K_o*Kd*a1*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 8*K_nt*Kd*O*a1*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*N*T*a1*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*N*O*T*a1*beta^2*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 16*K_nt*K_o*N*T*a1*beta^3*d^2 - 20*K_o*Kd*N*T*a1^2*beta^3*d - 16*K_nt*N*O*T*a1*beta^3*d^2 - 20*Kd*N*O*T*a1^2*beta^3*d - 8*K_o*Kd*T*a1*beta^2*d^2*delta + 8*N*O*T*a1*alpha_t*beta^2*d^2 - 8*Kd*O*T*a1*beta^2*d^2*delta + 8*K_o*N*T*a1*beta^2*d^2*m2 + 8*N*O*T*a1*beta^2*d^2*m2 - 20*K_o*Kd*N*T*a1^2*beta^2*d^2 - 16*K_o*Kd^2*N*T*a1^3*beta*d^2 - 16*K_o*Kd^2*N*T*a1^3*beta^2*d - 20*Kd*N*O*T*a1^2*beta^2*d^2 - 16*Kd^2*N*O*T*a1^3*beta*d^2 - 16*Kd^2*N*O*T*a1^3*beta^2*d - 4*K_o*Kd^2*T*a1^2*beta*d^2*delta - 4*K_o*Kd^2*T*a1^2*beta^2*d*delta - 8*K_o*N*T^2*a1*beta^2*d^2*delta + 8*K_nt*K_o*Kd*a1*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 4*Kd^2*O*T*a1^2*beta*d^2*delta - 4*Kd^2*O*T*a1^2*beta^2*d*delta - 8*N*O*T^2*a1*beta^2*d^2*delta + 8*K_nt*Kd*O*a1*beta^3*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*Kd*T*a1*beta*d^2*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*Kd*T*a1*beta^2*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*Kd*O*T*a1*beta*d^2*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*Kd*O*T*a1*beta^2*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*Kd*N*T*a1^2*beta*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*K_o*Kd*N*T*a1^2*beta^2*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*Kd*N*O*T*a1^2*beta*d^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 4*Kd*N*O*T*a1^2*beta^2*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2)
fO = O*alpha_o*beta^3 - O^2*beta^3*delta + K_o*beta^3*m3 + O*beta^3*m3 - O*alpha_o*beta^2*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - K_o*beta^2*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - O*beta^2*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - K_o*O*beta^3*delta + O*alpha_o*beta^2*d + K_o*beta^2*d*m3 + O*beta^2*d*m3 + O^2*beta^2*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - O^2*beta^2*d*delta - O*alpha_o*beta*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - K_o*beta*d*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - O*beta*d*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - Kd^2*O^2*a1^2*beta*delta - Kd^2*O^2*a1^2*d*delta + K_o*O*beta^2*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + O^2*beta*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 2*Kd*O*a1*alpha_o*beta^2 - K_o*O*beta^2*d*delta + 2*K_o*Kd*a1*beta^2*m3 + 4*K_nt*K_o*beta^2*d*m3 + 2*Kd*O*a1*beta^2*m3 + 4*K_nt*O*beta^2*d*m3 + Kd^2*O*a1^2*alpha_o*beta + Kd^2*O*a1^2*alpha_o*d - 2*Kd*O^2*a1*beta^2*delta - 4*K_nt*O^2*beta^2*d*delta + K_o*Kd^2*a1^2*beta*m3 + K_o*Kd^2*a1^2*d*m3 + Kd^2*O*a1^2*beta*m3 + Kd^2*O*a1^2*d*m3 - K_o*Kd^2*O*a1^2*beta*delta - K_o*Kd^2*O*a1^2*d*delta - Kd*O*a1*alpha_o*beta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - Kd*O*a1*alpha_o*d*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + K_o*O*beta*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - K_o*Kd*a1*beta*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - K_o*Kd*a1*d*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - Kd*O*a1*beta*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - Kd*O*a1*d*m3*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + 2*Kd*O*a1*alpha_o*beta*d + 2*K_o*Kd*a1*beta*d*m3 + 2*Kd*O*a1*beta*d*m3 + Kd*O^2*a1*beta*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + Kd*O^2*a1*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 2*K_o*Kd*O*a1*beta^2*delta - 4*K_nt*K_o*O*beta^2*d*delta - 2*Kd*O^2*a1*beta*d*delta - 2*N*O^2*T*a1*beta*d*delta + K_o*Kd*O*a1*beta*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) + K_o*Kd*O*a1*d*delta*(beta^3/(beta + d) + (beta^2*d)/(beta + d) + (2*Kd*a1*beta^2)/(beta + d) + (Kd^2*a1^2*beta)/(beta + d) + (Kd^2*a1^2*d)/(beta + d) + (2*Kd*a1*beta*d)/(beta + d) + (4*N*T*a1*beta*d)/(beta + d))^(1/2) - 2*K_o*Kd*O*a1*beta*d*delta + 2*N*O*T*a1*alpha_o*beta*d + 2*K_o*N*T*a1*beta*d*m3 + 2*N*O*T*a1*beta*d*m3 - 2*K_o*N*O*T*a1*beta*d*delta

nu = [fN, fT, fO]
@show nu = nu.subs(K_o , Kₒ )
@show nu = nu.subs(K_nt, Kₙₜ )
@show nu = nu.subs(Kd , Kd)
@show nu = nu.subs(a1,a₁)
@show nu = nu.subs(d ,d )
@show nu = nu.subs(a_nt ,aₙₜ)
@show nu = nu.subs(a_0  ,aₒ)
@show nu = nu.subs(alpha_t,αₜ)
@show nu = nu.subs(alpha_o,αₒ)
@show nu = nu.subs(alpha_n,αₙ)
@show nu = nu.subs(delta,δ)
@show nu = nu.subs(beta ,β)

display(nu)

nu





# --------- result below ----------
# ----- For matlab
# fN   = Dn1*alpha_n + Kd*NT*a1 - N*T*a1 - N*delta + m1
# fT   = Dt11*alpha_t + Kd*NT*a1 - N*T*a1 - T*delta + m2
# fx(3)   = -Dn0*O*a_0 + Dn1*K_o*a_0 - Do00*O*a_0 + Do01*K_o*a_0 - Do10*O*a_0 + Do11*K_o*a_0 + Do11*alpha_o - Dt00*O*a_0 + Dt01*K_o*a_0 - Dt10*O*a_0 + Dt11*K_o*a_0 - O*delta + m3
# fDo00 = -Do00*NT2*a_nt - Do00*O*a_0 + Do01*K_o*a_0 + Do10*K_nt*a_nt
# fDo10 = Do00*NT2*a_nt - Do10*K_nt*a_nt - Do10*O*a_0 + Do11*K_o*a_0
# fDo01 = Do00*O*a_0 - Do01*K_o*a_0 - Do01*NT2*a_nt + Do11*K_nt*a_nt
# fDo11 = Do01*NT2*a_nt + Do10*O*a_0 - Do11*K_nt*a_nt - Do11*K_o*a_0
# fDt00 = -Dt00*NT2*a_nt - Dt00*O*a_0 + Dt01*K_o*a_0 + Dt10*K_nt*a_nt
# fDt10 = Dt00*NT2*a_nt - Dt10*K_nt*a_nt - Dt10*O*a_0 + Dt11*K_o*a_0
# fDt01 = Dt00*O*a_0 - Dt01*K_o*a_0 - Dt01*NT2*a_nt + Dt11*K_nt*a_nt
# fDt11 = Dt01*NT2*a_nt + Dt10*O*a_0 - Dt11*K_nt*a_nt - Dt11*K_o*a_0
# fDn0  = -Dn0*O*a_0 + Dn1*K_o*a_0
# fDn1  = Dn0*O*a_0 - Dn1*K_o*a_0
# fNT  = -Kd*NT*a1 + N*T*a1 - NT^2*d - NT*beta + 2*NT2*d
# fNT2 = -Do00*NT2*a_nt - Do01*NT2*a_nt + Do10*K_nt*a_nt + Do11*K_nt*a_nt - Dt00*NT2*a_nt - Dt01*NT2*a_nt + Dt10*K_nt*a_nt + Dt11*K_nt*a_nt + 0.5*NT^2*d - NT2*beta - NT2*d








# ----- For mathematica
# fN   = Dn1*αₙ + Kd*NT*a₁ - N*T*a₁ - N*δ + m1
# fT   = Dt11*αₜ + Kd*NT*a₁ - N*T*a₁ - T*δ + m2
# fO   = -Dn0*O*aₒ + Dn1*Kₒ*aₒ - Do00*O*aₒ + Do01*Kₒ*aₒ - Do10*O*aₒ + Do11*Kₒ*aₒ + Do11*αₒ - Dt00*O*aₒ + Dt01*Kₒ*aₒ - Dt10*O*aₒ + Dt11*Kₒ*aₒ - O*δ + m3
# fDo00 = -Do00*NT2*aₙₜ - Do00*O*aₒ + Do01*Kₒ*aₒ + Do10*Kₙₜ*aₙₜ
# fDo10 = Do00*NT2*aₙₜ - Do10*Kₙₜ*aₙₜ - Do10*O*aₒ + Do11*Kₒ*aₒ
# fDo01 = Do00*O*aₒ - Do01*Kₒ*aₒ - Do01*NT2*aₙₜ + Do11*Kₙₜ*aₙₜ
# fDo11 = Do01*NT2*aₙₜ + Do10*O*aₒ - Do11*Kₒ*aₒ - Do11*Kₙₜ*aₙₜ
# fDt00 = -Dt00*NT2*aₙₜ - Dt00*O*aₒ + Dt01*Kₒ*aₒ + Dt10*Kₙₜ*aₙₜ
# fDt10 = Dt00*NT2*aₙₜ - Dt10*Kₙₜ*aₙₜ - Dt10*O*aₒ + Dt11*Kₒ*aₒ
# fDt01 = Dt00*O*aₒ - Dt01*Kₒ*aₒ - Dt01*NT2*aₙₜ + Dt11*Kₙₜ*aₙₜ
# fDt11 = Dt01*NT2*aₙₜ + Dt10*O*aₒ - Dt11*Kₒ*aₒ - Dt11*Kₙₜ*aₙₜ
# fDn0  = -Dn0*O*aₒ + Dn1*Kₒ*aₒ
# fDn1  = Dn0*O*aₒ - Dn1*Kₒ*aₒ
# fNT  = -Kd*NT*a₁ + N*T*a₁ - NT^2*d - NT*β + 2*NT2*d
# fNT2 = -Do00*NT2*aₙₜ - Do01*NT2*aₙₜ + Do10*Kₙₜ*aₙₜ + Do11*Kₙₜ*aₙₜ - Dt00*NT2*aₙₜ - Dt01*NT2*aₙₜ + Dt10*Kₙₜ*aₙₜ + Dt11*Kₙₜ*aₙₜ + 0.5*NT^2*d - NT2*d - NT2*β



























# ====================================================================================
# -------------------------- For demethylation model --------------------------------
using SymPy

@vars Kₒ Kd r₁ a₁ α₁ δₒ a0 β ζ₁ k₃ k₄ m1 m2 m3
@vars N  T  NT O  Dₒ Dᵒ D̄ Dₕ Dₜ Dᵗ Dₙ Dᴺ
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

@show exs = [dN, dT, dNT, dO, dDₒ, dDᵒ, dD̄, dDₕ, dDₜ, dDᵗ, dDₙ, dDᴺ]


# Convert variables
@vars Do DO D5mc Dh5mc Dt DT Dn DN
@show exs = exs.subs(Dₒ, Do)
@show exs = exs.subs(Dᵒ, DO)
@show exs = exs.subs(D̄, D5mc)
@show exs = exs.subs(Dₕ, Dh5mc)
@show exs = exs.subs(Dₜ, Dt)
@show exs = exs.subs(Dᵗ, DT)
@show exs = exs.subs(Dₙ, Dn)
@show exs = exs.subs(Dᴺ, DN)

# Convert Parameters

@vars K0 Kd r1 a1 alpha1 delta0 a0 beta xi1 k3 k4
@show exs = exs.subs(Kₒ, K0)
@show exs = exs.subs(Kd, Kd)
@show exs = exs.subs(r₁, r1)
@show exs = exs.subs(a₁, a1)
@show exs = exs.subs(α₁, alpha1)
@show exs = exs.subs(δₒ, delta0)
@show exs = exs.subs(a0, a0)
@show exs = exs.subs(β , beta)
@show exs = exs.subs(ζ₁, xi1)
@show exs = exs.subs(k₃, k3)
@show exs = exs.subs(k₄, k4)



# Result from the above convertion
dN = DN*K0*r1 + DN*alpha1 - Dn*N*r1 + Kd*NT*a1 - N*T*a1 - N*delta0 + m1
dT = DT*alpha1 + Kd*NT*a1 - N*T*a1 - T*delta0 + m2
dNT = -Kd*NT*a1 + N*T*a1
dO = DN*K0*xi1 + DO*K0*r1 + DO*alpha1 + DT*K0*r1 - Dn*O*xi1 - Do*O*r1 - Dt*O*r1 - O*delta0 + m3
dDo = DO*K0*r1 + Dh5mc*beta - Do*O*r1 - Do*a0
dDO = -DO*K0*r1 + Do*O*r1
dD5mc = -D5mc*NT*k3 + Do*a0
dDh5mc = D5mc*NT*k3 - Dh5mc*beta
dDt = DT*K0*r1 - Dt*NT*k4 - Dt*O*r1
dDT = -DT*K0*r1 + Dt*NT*k4 + Dt*O*r1
dDn = DN*K0*r1 + DN*K0*xi1 - Dn*N*r1 - Dn*O*xi1
dDN = -DN*K0*r1 - DN*K0*xi1 + Dn*N*r1 + Dn*O*xi1


# Matlab result N,T,O only after elimination
# With matlab variable
dN = (N*alpha1*r1 + O*alpha1*xi1 + K0*m1*r1 + N*m1*r1 + K0*m1*xi1 + O*m1*xi1 - N^2*delta0*r1 - K0*N*delta0*r1 - K0*N*delta0*xi1 - N*O*delta0*xi1 - K0*N*r1*xi1 + K0*O*r1*xi1)/(K0*r1 + N*r1 + K0*xi1 + O*xi1)
dT = m2 - T*delta0 + (alpha1*(Kd*O*r1 + N*T*k4))/(K0*Kd*r1 + Kd*O*r1 + N*T*k4)
dO = (Kd*(K0*a0*beta*m3 - K0*O*a0*beta*delta0 - K0*O*a0*beta*r1 - K0*O*a0*beta*xi1 + (K0^2*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (K0^2*a0*beta*r1*(Kd*O*r1 + N*T*k4))/(K0*Kd*r1 + Kd*O*r1 + N*T*k4) + (K0*O*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (N*O*T*alpha1*beta*k3)/Kd + (K0*N*T*a0*k3*m3)/Kd + (K0*N*T*beta*k3*m3)/Kd + (N*O*T*beta*k3*m3)/Kd + (K0*O*a0*beta*r1*(Kd*O*r1 + N*T*k4))/(K0*Kd*r1 + Kd*O*r1 + N*T*k4) - (N*O^2*T*beta*delta0*k3)/Kd - (N*O^2*T*beta*k3*r1)/Kd - (N*O^2*T*beta*k3*xi1)/Kd - (K0*N*O*T*a0*delta0*k3)/Kd - (K0*N*O*T*beta*delta0*k3)/Kd - (K0*N*O*T*a0*k3*r1)/Kd - (K0*N*O*T*beta*k3*r1)/Kd - (K0*N*O*T*a0*k3*xi1)/Kd - (K0*N*O*T*beta*k3*xi1)/Kd + (K0^2*N*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0^2*N*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (N*O^2*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0^2*N*T*a0*k3*r1*(Kd*O*r1 + N*T*k4))/(Kd*(K0*Kd*r1 + Kd*O*r1 + N*T*k4)) + (K0^2*N*T*beta*k3*r1*(Kd*O*r1 + N*T*k4))/(Kd*(K0*Kd*r1 + Kd*O*r1 + N*T*k4)) + (N*O^2*T*beta*k3*r1*(Kd*O*r1 + N*T*k4))/(Kd*(K0*Kd*r1 + Kd*O*r1 + N*T*k4)) + (K0*N*O*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (2*K0*N*O*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0*N*O*T*a0*k3*r1*(Kd*O*r1 + N*T*k4))/(Kd*(K0*Kd*r1 + Kd*O*r1 + N*T*k4)) + (2*K0*N*O*T*beta*k3*r1*(Kd*O*r1 + N*T*k4))/(Kd*(K0*Kd*r1 + Kd*O*r1 + N*T*k4))))/(K0*Kd*a0*beta + K0*N*T*a0*k3 + K0*N*T*beta*k3 + N*O*T*beta*k3)


# with julia variables
dN = (N*α₁*r₁ + O*α₁*ζ₁ + Kₒ*m1*r₁ + N*m1*r₁ + Kₒ*m1*ζ₁ + O*m1*ζ₁ - N^2*δₒ*r₁ - Kₒ*N*δₒ*r₁ - Kₒ*N*δₒ*ζ₁ - N*O*δₒ*ζ₁ - Kₒ*N*r₁*ζ₁ + Kₒ*O*r₁*ζ₁)/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)
dT = m2 - T*δₒ + (α₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)
dO = (Kd*(Kₒ*a0*β*m3 - Kₒ*O*a0*β*δₒ - Kₒ*O*a0*β*r₁ - Kₒ*O*a0*β*ζ₁ + (Kₒ^2*a0*β*ζ₁*(N*r₁ + O*ζ₁))/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁) + (Kₒ^2*a0*β*r₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄) + (Kₒ*O*a0*β*ζ₁*(N*r₁ + O*ζ₁))/(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁) + (N*O*T*α₁*β*k₃)/Kd + (Kₒ*N*T*a0*k₃*m3)/Kd + (Kₒ*N*T*β*k₃*m3)/Kd + (N*O*T*β*k₃*m3)/Kd + (Kₒ*O*a0*β*r₁*(Kd*O*r₁ + N*T*k₄))/(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄) - (N*O^2*T*β*δₒ*k₃)/Kd - (N*O^2*T*β*k₃*r₁)/Kd - (N*O^2*T*β*k₃*ζ₁)/Kd - (Kₒ*N*O*T*a0*δₒ*k₃)/Kd - (Kₒ*N*O*T*β*δₒ*k₃)/Kd - (Kₒ*N*O*T*a0*k₃*r₁)/Kd - (Kₒ*N*O*T*β*k₃*r₁)/Kd - (Kₒ*N*O*T*a0*k₃*ζ₁)/Kd - (Kₒ*N*O*T*β*k₃*ζ₁)/Kd + (Kₒ^2*N*T*a0*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ^2*N*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (N*O^2*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ^2*N*T*a0*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (Kₒ^2*N*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (N*O^2*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (Kₒ*N*O*T*a0*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (2*Kₒ*N*O*T*β*k₃*ζ₁*(N*r₁ + O*ζ₁))/(Kd*(Kₒ*r₁ + N*r₁ + Kₒ*ζ₁ + O*ζ₁)) + (Kₒ*N*O*T*a0*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄)) + (2*Kₒ*N*O*T*β*k₃*r₁*(Kd*O*r₁ + N*T*k₄))/(Kd*(Kₒ*Kd*r₁ + Kd*O*r₁ + N*T*k₄))))/(Kₒ*Kd*a0*β + Kₒ*N*T*a0*k₃ + Kₒ*N*T*β*k₃ + N*O*T*β*k₃)
