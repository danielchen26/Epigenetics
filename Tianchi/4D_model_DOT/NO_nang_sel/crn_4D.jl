# ================ CRN ======================
using DiffEqBiological, LinearAlgebra
using Plots;gr()
using DataFrames, Queryverse, Latexify



#  ----------      🔺        +          💠      CRN model ------------
#  -------- Demethylation   +   TF 2sites binding -------------------

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

end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn kh beta m1 m2 m3


p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]

@add_constraints Demethy_TF_MC begin
  Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end

ss = steady_states(Demethy_TF_MC,p)
stability(ss,Demethy_TF_MC,p)

bif_grid_dia = bifurcation_grid_diagram(Demethy_TF_MC, p, :kh, 0.:0.1:0.55,  :beta, (0.1,1.))
bp = plot(bif_grid_dia)

params = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
ss = steady_states(Demethy_TF_MC,params)
stability(ss,Demethy_TF_MC,params)

dfc = DataFrame(vcat(ss))
dfc.name = Demethy_TF_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame
