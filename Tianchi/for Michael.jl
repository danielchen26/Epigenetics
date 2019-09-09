using SymPy

@vars N T NT NT2 D00_O D10_O D01_O D11_O O D00_T D10_T D01_T D11_T D0_N D1_N K_o K_nt Kd a1 d a_nt a0 alpha_t alpha_o alpha_n delta beta m1 m2 m3
# ------ Symstem of equations
fN   = -N*delta + alpha_n*D1_N - N*T*a1 + NT*Kd*a1 + m1
fT   = -T*delta + alpha_t*D11_T - N*T*a1 + NT*Kd*a1 + m2
fO   = -O*delta + alpha_o*D11_O + K_o*a0*D01_O + K_o*a0*D01_T + K_o*a0*D1_N + K_o*a0*D11_O + K_o*a0*D11_T - O*a0*D0_N - O*a0*D00_O - O*a0*D00_T - O*a0*D10_O - O*a0*D10_T + m3
fD00_O = K_o*a0*D01_O - NT2*a_nt*D00_O - O*a0*D00_O + a_nt*K_nt*D10_O
fD10_O = K_o*a0*D11_O + NT2*a_nt*D00_O - O*a0*D10_O - a_nt*K_nt*D10_O
fD01_O = -K_o*a0*D01_O - NT2*a_nt*D01_O + O*a0*D00_O + a_nt*K_nt*D11_O
fD11_O = -K_o*a0*D11_O + NT2*a_nt*D01_O + O*a0*D10_O - a_nt*K_nt*D11_O
fD00_T = K_o*a0*D01_T - NT2*a_nt*D00_T - O*a0*D00_T + a_nt*K_nt*D10_T
fD10_T = K_o*a0*D11_T + NT2*a_nt*D00_T - O*a0*D10_T - a_nt*K_nt*D10_T
fD01_T = -K_o*a0*D01_T - NT2*a_nt*D01_T + O*a0*D00_T + a_nt*K_nt*D11_T
fD11_T = -K_o*a0*D11_T + NT2*a_nt*D01_T + O*a0*D10_T - a_nt*K_nt*D11_T
fD0_N  = K_o*a0*D1_N - O*a0*D0_N
fD1_N  = -K_o*a0*D1_N + O*a0*D0_N
fNT  = -NT*beta - d*NT^2 + 2*d*NT2 + N*T*a1 - NT*Kd*a1
fNT2 = (1/2)*d*NT^2 - d*NT2 - beta*NT2 - NT2*a_nt*D00_O - NT2*a_nt*D00_T - NT2*a_nt*D01_O - NT2*a_nt*D01_T + a_nt*K_nt*D10_O + a_nt*K_nt*D10_T + a_nt*K_nt*D11_O + a_nt*K_nt*D11_T
@show exs = [ fN, fT, fO, fD00_O, fD10_O, fD01_O, fD11_O, fD00_T, fD10_T, fD01_T, fD11_T, fD0_N, fD1_N, fNT, fNT2]



# 1. Substitute with conservation law
D11_O_sub = 1 - D00_O - D10_O - D01_O
D11_T_sub = 1 - D00_T - D10_T - D01_T
D1_N_sub  = 1 - D0_N
@show exs = exs.subs.(D1_N, D1_N_sub )
@show exs =exs.subs.(D11_T, D11_T_sub )
@show exs =exs.subs.(D11_O, D11_O_sub )
@. exs = simplify(exs)

# 2. Solve for other variables except for 1:=> N, 2:=>T, 3:=>O.
dD00_O = solve(fD00_O, D00_O)
@show exs =exs.subs.(D00_O, dD00_O)[1]
dD10_O = solve(fD10_O, D10_O)
@show exs =exs.subs.(D10_O, dD10_O)[1]
dD01_O = solve(fD01_O, D01_O)
@show exs =exs.subs.(D01_O, dD01_O)[1]
@. exs = simplify(exs)

dD00_T = solve(fD00_T, D00_T)
@show exs =exs.subs.(D00_T, dD00_T)[1]
dD10_T = solve(fD10_T, D10_T)
@show exs =exs.subs.(D10_T, dD10_T)[1]
dD01_T = solve(fD01_T, D01_T)
@show exs =exs.subs.(D01_T, dD01_T)[1]
@. exs = simplify(exs)

dD0_N = solve(fD0_N, D0_N)
@show exs =exs.subs.(D0_N, dD0_N)[1]
@. exs = simplify(exs)

dNT = solve(fNT, NT)
@show exs =exs.subs.(NT, dNT)[1]
dNT2 = solve(fNT2, NT2)
@show exs =exs.subs.(NT2, dNT2)[1]
@. exs = simplify(exs)

@show exs
