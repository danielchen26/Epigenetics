using DiffEqBiological, LinearAlgebra, BlackBoxOptim, StaticArrays
using DataFrames, Queryverse, Suppressor
using Plots

# ========= Paper model ------- (matlab syntax: C)
ps_MC = @reaction_network begin
    # N T complex
    (a1, Kd*a1),               N + T ↔ NT
    (d, d),                   NT + NT ↔ NT2

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT2 + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT2 + Do01 ↔ Do11
    (ao, Ko*ao),              O + Do00 ↔ Do01
    (ao, Ko*ao),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT2 + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT2 + Dt01 ↔ Dt11
    (ao, Ko*ao),              O + Dt00 ↔ Dt01
    (ao, Ko*ao),              O + Dt10 ↔ Dt11
    # Nanog
    (ao, Ko*ao),              O + Dn0 ↔ Dn1

    # --- Protein production
    alphaT,                  Dt11 → Dt11 + T
    alphaO,                  Do11 → Do11 + O
    alphaN,                  Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),    (N, T, O) → ∅
    #(β,β),             (NT, NT2) → ∅

    # ---- NTO rate control m1 m2 m3-----
    m1,                       ∅ → N
    m2,                       ∅ → T
    m3,                       ∅ → O
end Ko K_nt Kd a1 d a_nt ao alphaT alphaO alphaN delta m1 m2 m3





#   ------ Demethylation CRN ---------
Demethy_crn_MC = @reaction_network begin
    # N T complex
    (a1, Kd*a1),       N + T ↔ NT              # T is recruited by N

    # --- Promoter binding
    # Oct4 cycle
    (r1, K0*r1 ),      O + Do ↔ DO              # O Auto-activation
    alpha1,            DO → DO + O              # O translation
    delta0,            O → ∅                    # Deg
    # Oct4 de-Methylation cycle
    a0,                Do → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    beta,              D5hmc →  Do              # 5hmC -> C by dilution(replication)
    # TET protein
    (r1,K0*r1),        Dt + O ↔ DT              # O -> T promoter
    alpha1,            DT → DT + T              # T translation
    delta0,            T → ∅                    # Deg
    # Nanog cycle
    (r1, K0*r1 ),       N + Dn ↔ DN             # N Auto-activation
    alpha1,             DN → DN + N             # N translation
    delta0,             N → ∅                   # Deg
    # O activate N
    (xi1, K0*xi1),      O + Dn ↔ DN             # O -> N promoter

    # ---- Nanog guided Tet1/2
    k3,                NT + D5mc → D5hmc + NT   # NT oxidize 5mc -> 5hmC
    # k₄,               NT + Dt → DT + NT       # NT induce not only Oct4 but T1

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end K0 Kd r1 a1 alpha1 delta0 a0 beta xi1 k3 m1 m2 m3

@add_constraints Demethy_crn_MC begin
  Do + DO + D5mc + D5hmc  = 1
  Dt + DT = 1
  Dn + DN = 1
end



#  -------- Demethylation   +   TF 2sites binding ------------------
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    (d, d),                   NT + NT ↔ NT2     #💚 Need dimerization?

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT2 + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT2 + Do01 ↔ Do11
    (aO, KO*aO),              O + Do00 ↔ Do01
    (aO, KO*aO),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT2 + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT2 + Dt01 ↔ Dt11
    (aO, KO*aO),              O + Dt00 ↔ Dt01
    (aO, KO*aO),              O + Dt10 ↔ Dt11
    # Nanog
    (aO, KO*aO),              O + Dn0 ↔ Dn1

    # --- Protein production
    alphaT,                   Dt11 → Dt11 + T
    alphaO,                   Do11 → Do11 + O
    alphaN,                   Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),           (N, T, O) → ∅

    # ============   🔺 ===============
    # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    kh,                  NT + D5mc → D5hmc + NT     # NT oxidize 5mc -> 5hmC
    beta,                D5hmc →  Do00              # 5hmC -> C by

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn beta kh m1 m2 m3
    # ============   🔺 ===============  Demethy_crn_MatlabC
    # --- Promoter binding
    # # Oct4 cycle
    # (r1, K0*r1 ),      O + Do ↔ DO              # O Auto-activation
    # alpha1,            DO → DO + O              # O translation
    # delta0,            O → ∅                    # Deg
    #
    # # TET protein
    # (r1,K0*r1),        Dt + O ↔ DT              # O -> T promoter
    # alpha1,            DT → DT + T              # T translation
    # delta0,            T → ∅                    # Deg
    #
    # # Nanog cycle ----💚 need N auto-activation?
    # (r1, K0*r1 ),       N + Dn ↔ DN             # N Auto-activation
    # alpha1,             DN → DN + N             # N translation
    # delta0,             N → ∅                   # Deg
    # O activate N
    # (xi1, K0*xi1),      O + Dn ↔ DN             # O -> N promoter
    #
    # ---- Nanog guided Tet1/2
    #
    # k₄,               NT + Dt → DT + NT       # NT induce not only Oct4 but T1
p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
@add_constraints Demethy_TF_MC begin
  Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end







# --------------- Global optimizer for same SSS in paper model--------------
# for 🔺 model -----
params = [0.3, 0.1, 1, 1, 1, 1, 0.5, 1, 1, 1,    0., 0.05, 0.]
ss =  steady_states(Demethy_crn_MC,params)

function rn_SSS_match(p)
    params[3:10] = p
    # @show params
    ss = steady_states(Demethy_crn_MC,params)
    target = [0.6954, 0.7349,  0.6849]
    NTO_idx = [1,2,4]
    # @show length(ss)
    cost = [norm(ss[i][[1,2,4]] .- target) for i = 1: length(ss)]
    # @show cost
  return minimum(cost)
end
@time res = bboptimize(rn_SSS_match; SearchRange = (0.0, 3.0), NumDimensions = 8)





# --------------- Global optimizer for same SSS in Demethy_TF model--------------
# KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn beta kh m1 m2 m3
params = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
ss =  steady_states(Demethy_TF_MC,params)

function rn_SSS_match2(p)
    # params[[4:5; 8:14]] = p
    params[12:14] = p
    # @show params
    @suppress ss = steady_states(Demethy_TF_MC,params)
    target = [0.6954, 0.7349,  0.6849]
    NTO_idx = [1,2,4]
    # @show length(ss)
    cost = [norm(ss[i][[1,2,4]] .- target) for i = 1: length(ss)]
    # @show cost
    cc =  isempty(cost) == false ? minimum(cost) : 1e3
    return cc
end



rg1 = (0.,50.)#; rg2= (100., 10000.)
# range = [rg1,rg1,rg2,rg2,rg1,rg1,rg1,rg1,rg1,rg1,rg1]
@time res = bboptimize(rn_SSS_match2; SearchRange = rg1, NumDimensions =3)
@time res = compare_optimizers(rn_SSS_match2; SearchRange = rg1, NumDimensions =9)

best_fitness(res)






params =[0.3, 0.2, 0.1, 8.714, 7.42169, 1781.88, 1208.44, 9.91857, 1.66066, 3.54602, 1.89894, 9.15457, 1.90002, 9.96402, 0., 0.05, 0.]
ss =  steady_states(Demethy_TF_MC,params)
stability(ss, Demethy_TF_MC,params)
dfc = DataFrame(vcat(ss))
dfc.name = Demethy_TF_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame









# --------- find the term diff between ps_MC & Demethy_TF_MC --------
# ps_MC ode
f1_N    = :(m1 + Dn1*alphaN - N*delta - N*T*a1 + NT*Kd*a1)
f1_T    = :(m2 + Dt11*alphaT - T*delta - N*T*a1 + NT*Kd*a1)
f1_NT   = :(-d*NT^2 + 2*d*NT2 + N*T*a1 - NT*Kd*a1)
f1_NT2  = :((1/2)*d*NT^2 - d*NT2 + K_nt*a_nt*Do10 + K_nt*a_nt*Do11 + K_nt*a_nt*Dt10 + K_nt*a_nt*Dt11 - NT2*a_nt*Do00 - NT2*a_nt*Do01 - NT2*a_nt*Dt00 - NT2*a_nt*Dt01)
f1_Do00 = :(K_nt*a_nt*Do10 + Ko*ao*Do01 - NT2*a_nt*Do00 - O*ao*Do00)
f1_Do10 = :(-K_nt*a_nt*Do10 + Ko*ao*Do11 + NT2*a_nt*Do00 - O*ao*Do10)
f1_Do01 = :(K_nt*a_nt*Do11 - Ko*ao*Do01 - NT2*a_nt*Do01 + O*ao*Do00)
f1_Do11 = :(-K_nt*a_nt*Do11 - Ko*ao*Do11 + NT2*a_nt*Do01 + O*ao*Do10)
f1_O    = :(m3 + Do11*alphaO - O*delta + Ko*ao*Dn1 + Ko*ao*Do01 + Ko*ao*Do11 + Ko*ao*Dt01 + Ko*ao*Dt11 - O*ao*Dn0 - O*ao*Do00 - O*ao*Do10 - O*ao*Dt00 - O*ao*Dt10)
f1_Dt00 = :(K_nt*a_nt*Dt10 + Ko*ao*Dt01 - NT2*a_nt*Dt00 - O*ao*Dt00)
f1_Dt10 = :(-K_nt*a_nt*Dt10 + Ko*ao*Dt11 + NT2*a_nt*Dt00 - O*ao*Dt10)
f1_Dt01 = :(K_nt*a_nt*Dt11 - Ko*ao*Dt01 - NT2*a_nt*Dt01 + O*ao*Dt00)
f1_Dt11 = :(-K_nt*a_nt*Dt11 - Ko*ao*Dt11 + NT2*a_nt*Dt01 + O*ao*Dt10)
f1_Dn0  = :(Ko*ao*Dn1 - O*ao*Dn0)
f1_Dn1  = :(-Ko*ao*Dn1 + O*ao*Dn0)


# Demethy_TF_MC ode
f2_N     = :(m1 + Dn1*alphaN - N*delta - N*T*a1 + NT*Kd*a1)
f2_T     = :(m2 + Dt11*alphaT - T*delta - N*T*a1 + NT*Kd*a1)
f2_NT    = :(-d*NT^2 + 2*d*NT2 + N*T*a1 - NT*Kd*a1)
f2_NT2   = :((1/2)*d*NT^2 - d*NT2 + K_nt*a_nt*Do10 + K_nt*a_nt*Do11 + K_nt*a_nt*Dt10 + K_nt*a_nt*Dt11 - NT2*a_nt*Do00 - NT2*a_nt*Do01 - NT2*a_nt*Dt00 - NT2*a_nt*Dt01)
f2_Do00  = :(-a_dn*Do00 + beta*D5hmc + KO*aO*Do01 + K_nt*a_nt*Do10 - NT2*a_nt*Do00 - O*aO*Do00)
f2_Do10  = :(KO*aO*Do11 - K_nt*a_nt*Do10 + NT2*a_nt*Do00 - O*aO*Do10)
f2_Do01  = :(K_nt*a_nt*Do11 - Ko*ao*Do01 - NT2*a_nt*Do01 + O*ao*Do00)
f2_Do11  = :(-KO*aO*Do11 - K_nt*a_nt*Do11 + NT2*a_nt*Do01 + O*aO*Do10)
f2_O     = :(m3 + Do11*alphaO - O*delta + KO*aO*Dn1 + KO*aO*Do01 + KO*aO*Do11 + KO*aO*Dt01 + KO*aO*Dt11 - O*aO*Dn0 - O*aO*Do00 - O*aO*Do10 - O*aO*Dt00 - O*aO*Dt10)
f2_Dt00  = :(KO*aO*Dt01 + K_nt*a_nt*Dt10 - NT2*a_nt*Dt00 - O*aO*Dt00)
f2_Dt10  = :(KO*aO*Dt11 - K_nt*a_nt*Dt10 + NT2*a_nt*Dt00 - O*aO*Dt10)
f2_Dt01  = :(-KO*aO*Dt01 + K_nt*a_nt*Dt11 - NT2*a_nt*Dt01 + O*aO*Dt00)
f2_Dt11  = :(-KO*aO*Dt11 - K_nt*a_nt*Dt11 + NT2*a_nt*Dt01 + O*aO*Dt10)
f2_Dn0   = :(KO*aO*Dn1 - O*aO*Dn0)
f2_Dn1   = :(-KO*aO*Dn1 + O*aO*Dn0)
f2_D5mc  = :(a_dn*Do00 - kh*NT*D5mc)
f2_D5hmc = :(-beta*D5hmc + kh*NT*D5mc)


# Algebra.:-(f1_N,f2_N)
@force using Reduce.Algebra
fd_N     =   f1_N     -   f2_N
fd_T     =   f1_T     -   f2_T
fd_NT    =   f1_NT    -   f2_NT
fd_NT2   =   f1_NT2   -   f2_NT2
fd_Do00  =   f1_Do00  -   f2_Do00
fd_Do10  =   f1_Do10  -   f2_Do10
fd_Do01  =   f1_Do01  -   f2_Do01
fd_Do11  =   f1_Do11  -   f2_Do11
fd_O     =   f1_O     -   f2_O
fd_Dt00  =   f1_Dt00  -   f2_Dt00
fd_Dt10  =   f1_Dt10  -   f2_Dt10
fd_Dt01  =   f1_Dt01  -   f2_Dt01
fd_Dt11  =   f1_Dt11  -   f2_Dt11
fd_Dn0   =   f1_Dn0   -   f2_Dn0
fd_Dn1   =   f1_Dn1   -   f2_Dn1
f2_D5mc
f2_D5hmc
# :(f1_N = $f1_N) |> rcall






@vars x y z
