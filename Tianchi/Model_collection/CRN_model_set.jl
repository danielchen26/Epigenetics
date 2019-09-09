using DiffEqBiological, DifferentialEquations
using DataFrames, LinearAlgebra, Plots
using Queryverse

# # ========= Paper model ------- (matlab syntax: NC)
# ps = @reaction_network begin
#     # N T complex
#     (a₁, Kd*a₁),         N + T ↔ NT
#     (d, d),             NT + NT ↔ NT2
#
#     # --- Promoter binding
#     # Oct4
#     (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₀ᴼ ↔ D₁₀ᴼ
#     (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₁ᴼ ↔ D₁₁ᴼ
#     (aₒ, Kₒ*aₒ),          O + D₀₀ᴼ ↔ D₀₁ᴼ
#     (aₒ, Kₒ*aₒ),          O + D₁₀ᴼ ↔ D₁₁ᴼ
#     # TET
#     (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₀ᵀ ↔ D₁₀ᵀ
#     (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₁ᵀ ↔ D₁₁ᵀ
#     (aₒ, Kₒ*aₒ),          O + D₀₀ᵀ ↔ D₀₁ᵀ
#     (aₒ, Kₒ*aₒ),          O + D₁₀ᵀ ↔ D₁₁ᵀ
#     # Nanog
#     (aₒ, Kₒ*aₒ),          O + D₀ᴺ ↔ D₁ᴺ
#
#     # --- Protein production
#     αₜ,                 D₁₁ᵀ → D₁₁ᵀ + T
#     αₒ,                 D₁₁ᴼ → D₁₁ᴼ + O
#     αₙ,                 D₁ᴺ → D₁ᴺ + N
#
#     # --- Dilution and Degradation
#     (δ,δ,δ),           (N, T, O) → ∅
#     #(β,β),             (NT, NT2) → ∅
# end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β
#
# p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0.]
# @add_constraints ps begin
#   D₀₀ᴼ + D₁₀ᴼ + D₀₁ᴼ + D₁₁ᴼ = 1
#   D₀₀ᵀ + D₁₀ᵀ + D₀₁ᵀ + D₁₁ᵀ = 1
#   D₀ᴺ  + D₁ᴺ = 1
# end
# ss2 = steady_states(ps,p)
# stability(ss2[2],ps,p)


# # ========== Methylation CRN ========= (matlab syntax: NC)
# Demethy_crn = @reaction_network begin
#     # N T complex
#     (a₁, Kd*a₁),       N + T ↔ NT        # T is recruited by N
#
#     # --- Promoter binding
#     # Oct4 cycle
#     (r₁, Kₒ*r₁ ),      O + Dₒ ↔ Dᵒ           # O Auto-activation
#     α₁,               Dᵒ → Dᵒ + O           # O translation
#     δₒ,               O → ∅                 # Deg
#     # Oct4 de-Methylation cycle
#     a0,               Dₒ → D̄               # D + DNMT ↔ C₁ → D̄ + DNMT
#     β,                Dₕ →  Dₒ               # 5hmC -> C by dilution(replication)
#     # TET protein
#     (r₁,Kₒ*r₁),       Dₜ + O ↔ Dᵗ            # O -> T promoter
#     α₁,               Dᵗ → Dᵗ + T            # T translation
#     δₒ,               T → ∅                 # Deg
#     # Nanog cycle
#     (r₁, Kₒ*r₁ ),       N + Dₙ ↔ Dᴺ            # N Auto-activation
#     α₁,                Dᴺ → Dᴺ + N            # N translation
#     δₒ,                N → ∅                  # Deg
#     # O activate N
#     (ζ₁, Kₒ*ζ₁),       O + Dₙ ↔ Dᴺ            # O -> N promoter
#
#     # ---- Nanog guided Tet1/2
#     k₃,               NT + D̄ → Dₕ + NT       # NT oxidize 5mc -> 5hmC
#     # k₄,               NT + Dₜ → Dᵗ + NT       # NT induce not only Oct4 but T1
# end Kₒ Kd r₁ a₁ α₁ δₒ a0 β ζ₁ k₃# k₄
# @add_constraints Demethy_crn begin
#   Dₒ + Dᵒ + D̄ + Dₕ  = 1
#   Dₜ + Dᵗ = 1
#   Dₙ + Dᴺ = 1
# end
#
# params = [0.3, 0.1, 1, 1, 1, 1, 0.5, 1, 1, 1]
# ss = steady_states(Demethy_crn,params)
# stability(ss,Demethy_crn,params)
#
# df1 = DataFrame(vcat(ss))
# df1.name = Demethy_crn.syms
# var = [:N, :T, :O]
# @show df11 = df1 |> @filter(_.name in var) |> DataFrame
#
#
# bif_grid_dia = bifurcation_grid_diagram(Demethy_crn, params, :k₄, 1.:1.:10., :k₃, (1e-2, 5))
# plot(bif_grid_dia)



# ================= Below syntax is Maltab Compatible  ================================
# ================= Below syntax is Maltab Compatible  ================================



# ========= Paper model ------- (matlab syntax: C)
ps_MC = @reaction_network begin
    # N T complex
    (a1, Kd*a1),         N + T ↔ NT
    (d, d),             NT + NT ↔ NT2

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT2 + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT2 + Do01 ↔ Do11
    (ao, Ko*ao),          O + Do00 ↔ Do01
    (ao, Ko*ao),          O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT2 + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT2 + Dt01 ↔ Dt11
    (ao, Ko*ao),          O + Dt00 ↔ Dt01
    (ao, Ko*ao),          O + Dt10 ↔ Dt11
    # Nanog
    (ao, Ko*ao),          O + Dn0 ↔ Dn1

    # --- Protein production
    alphaT,                 Dt11 → Dt11 + T
    alphaO,                 Do11 → Do11 + O
    alphaN,                 Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),           (N, T, O) → ∅
    #(β,β),             (NT, NT2) → ∅

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O

end Ko K_nt Kd a1 d a_nt ao alphaT alphaO alphaN delta m1 m2 m3

p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1,    0., 0.05, 0.]
@add_constraints ps_MC begin
  Do00 + Do10 + Do01 + Do11 = 1
  Dt00 + Dt10 + Dt01 + Dt11 = 1
  Dn0  + Dn1 = 1
end
ss = steady_states(ps_MC,p)
@show stability(ss,ps_MC,p)


dfc = DataFrame(vcat(ss))
dfc.name = ps_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame


ps_MC.f_symfuncs

# fN   = Dn1*alpha_n + Kd*NT*a1 - N*T*a1 - N*delta + m1
# fT   = Dt11*alpha_t + Kd*NT*a1 - N*T*a1 - T*delta + m2
# fO   = -Dn0*O*a_0 + Dn1*K_o*a_0 - Do00*O*a_0 + Do01*K_o*a_0 - Do10*O*a_0 + Do11*K_o*a_0 + Do11*alpha_o - Dt00*O*a_0 + Dt01*K_o*a_0 - Dt10*O*a_0 + Dt11*K_o*a_0 - O*delta + m3
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
# ConsDo = 1 - Do10 - Do01 - Do00 - Do11;
# ConsDt = 1 - Dt10 - Dt01 - Dt00 - Dt11;
# ConsDn = 1 - Dn0 - Dn1;








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
end K0 Kd r1 a1 alpha1 delta0 a0 beta xi1 k3
