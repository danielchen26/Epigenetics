# ================ CRN ======================
using DiffEqBiological, LinearAlgebra
using Plots;gr()
using DataFrames, Queryverse, Latexify
include(pwd()*"/functions.jl")



# ======= Original CRN ========
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
    # (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
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
@add_constraints Demethy_TF_MC begin
  Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1 end

p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
p = [0.1, 0.2, 0.2, 5., 0., 500., 500., 5.0, 5.0, 5.0, 5., 5., 5., 5.,  0., 0.05, 0.]


ss = steady_states(Demethy_TF_MC,p)
sb = stability(ss,Demethy_TF_MC,p)
sb2 = stability_tianchi(ss,Demethy_TF_MC,p, 3)

@show round_ss = [round.(i,digits = 2) for i in ss]
@show [round.(real(JE_stability(i, Demethy_TF_MC, p)[2]),digits = 2) for i in ss]





plotly()
# Bifurcation diagram
p = [0.24, 0.19,0.26,2.7,5.0,500.,500.,5.0,5.0,5.0,5.0,5.0,5.0,5.0,  0., 0.05, 0.]
bif_grid_dia = bifurcation_grid_diagram(Demethy_TF_MC, p, :KO, 0.:0.1:0.55,  :Kd, (0.1,.3))
bp = plot(bif_grid_dia)
plot!(bp,[[],[],[],[]],color=[:blue :cyan :orange :red],label = ["Stable Real" "Stable Complex" "Unstable Complex" "Unstable Real"])

bif = bifurcation_grid_2d(Demethy_TF_MC, p, :KO, 0.:0.1:0.55,  :Kd, 0.1:0.1:1.)
plot(bif)
scatter!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])



dfc = DataFrame(vcat(ss))
dfc.name = Demethy_TF_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame






# First Visulization =====================
# ===== 3d model DOT defined by N-T-O-----

using Interact #, Suppressor, ProgressMeter
function Interactive_model(Demethy_TF_MC)
    @manipulate for KO = 0:0.01:1.0, K_nt = 0:0.01:1.0, Kd = 0:0.01:1.0, a1 = 0:0.1:10.0, d = 0:0.1:10.0, a_nt = 0:10:1000.0, aO = 0:10:1000.0, alphaT = 0:0.1:10.0, alphaO = 0:0.1:10.0, alphaN = 0:0.1:10.0, delta = 0:0.1:10.0, a_dn = 0:0.1:10.0, kh = 0:0.1:10.0, beta= 0:0.1:10.0
        p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta, 0., 0.05, 0.]
        ss = steady_states(Demethy_TF_MC,p)
        sort!(ss, by = x -> x[1])
        sb = stability(ss,Demethy_TF_MC,p)

        dfc = DataFrame(vcat(ss))
        dfc.name = Demethy_TF_MC.syms
        var = [:N, :T, :O]
        df_iPS = dfc |> @filter(_.name in var) |> DataFrame
        @show Matrix(df_iPS)

        # Need to modify below
        # if length(ss) >2
        #     DOT   = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3]))
        #     @show DOT
        #
        # else
        #      @show "single state"
        # end

        # eigs = [round.(JE_stability(i,Demethy_TF_MC,p)[2] ,digits = 2) for i in ss]
        # println("The eigs are:")
        # @show eigs
        # println("Eig values:",vcat(eigs...))

        eigs = [round.(real(JE_stability(i, Demethy_TF_MC, p)[2]),digits = 2) for i in ss]
        @show eigs # show eigen values real part 2digits
        @show string.(sb)
        plot(sort([i[[1,2,9]] for i in ss]), xlabel = "N-T-O", label =string.(sb), marker = (:hexagon, 10, 0.7, stroke(1, 0.1, :black, :dot)))
        # println("The SS states are:")
        # sort!(ss, by = x -> x[1])
        # @show ss
        # @show [round.(EigenVector(Demethy_TF_MC,p,i),digits = 2) for i in ss]
    end
end

# ==========================================================================
# =================== Different_Toplogies ==================================

# 1. NT2 regulate 💠, NT regulate 🔺, Nanog self: NONE {Original CRN Model}
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
    # (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
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

# 2. NT regulate 💠, NT2 regulate 🔺, Nanog self: NONE
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    (d, d),                   NT + NT ↔ NT2     #💚 Need dimerization?

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT + Do01 ↔ Do11
    (aO, KO*aO),              O + Do00 ↔ Do01
    (aO, KO*aO),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT + Dt01 ↔ Dt11
    (aO, KO*aO),              O + Dt00 ↔ Dt01
    (aO, KO*aO),              O + Dt10 ↔ Dt11
    # Nanog
    (aO, KO*aO),              O + Dn0 ↔ Dn1
    # (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
    # --- Protein production
    alphaT,                   Dt11 → Dt11 + T
    alphaO,                   Do11 → Do11 + O
    alphaN,                   Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),           (N, T, O) → ∅

    # ============   🔺 ===============
    # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    kh,                  NT2 + D5mc → D5hmc + NT2     # NT oxidize 5mc -> 5hmC
    beta,                D5hmc →  Do00              # 5hmC -> C by

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn kh beta m1 m2 m3

# 3. NT2 regulate 💠, NT2 regulate 🔺, Nanog self: NONE
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
    # (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
    # --- Protein production
    alphaT,                   Dt11 → Dt11 + T
    alphaO,                   Do11 → Do11 + O
    alphaN,                   Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),           (N, T, O) → ∅

    # ============   🔺 ===============
    # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    kh,                  NT2 + D5mc → D5hmc + NT2     # NT oxidize 5mc -> 5hmC
    beta,                D5hmc →  Do00              # 5hmC -> C by

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn kh beta m1 m2 m3

# 4. NT regulate 💠,   NT regulate 🔺, Nanog self: NONE {rate d = 0}
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    (d, d),                   NT + NT ↔ NT2     #💚 Need dimerization?

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT + Do01 ↔ Do11
    (aO, KO*aO),              O + Do00 ↔ Do01
    (aO, KO*aO),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT + Dt01 ↔ Dt11
    (aO, KO*aO),              O + Dt00 ↔ Dt01
    (aO, KO*aO),              O + Dt10 ↔ Dt11
    # Nanog
    (aO, KO*aO),              O + Dn0 ↔ Dn1
    # (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
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

# 5. NT regulate 💠,   NT regulate 🔺, Nanog self: Yes {rate d = 0}
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    (d, d),                   NT + NT ↔ NT2     #💚 Need dimerization?

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT + Do01 ↔ Do11
    (aO, KO*aO),              O + Do00 ↔ Do01
    (aO, KO*aO),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT + Dt01 ↔ Dt11
    (aO, KO*aO),              O + Dt00 ↔ Dt01
    (aO, KO*aO),              O + Dt10 ↔ Dt11
    # Nanog
    (aO, KO*aO),              O + Dn0 ↔ Dn1
    (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
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

# 6. NT regulate 💠,   NT regulate 🔺, Nanog Dimer self: Yes {d is for N2}
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    (d, d),                   N + N ↔ N2     #💚 Need dimerization?

    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT + Do01 ↔ Do11
    (aO, KO*aO),              O + Do00 ↔ Do01
    (aO, KO*aO),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT + Dt01 ↔ Dt11
    (aO, KO*aO),              O + Dt00 ↔ Dt01
    (aO, KO*aO),              O + Dt10 ↔ Dt11
    # Nanog
    (aO, KO*aO),              O + Dn0 ↔ Dn1
    (aO, KO*aO),              N2 + Dn0 ↔ Dn1 #  a self regulation is added
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

# 7. NT2 regulate 💠,   NT regulate 🔺, Nanog  self: Yes
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
    (aO, KO*aO),              N + Dn0 ↔ Dn1 #  a self regulation is added
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


# ==== Add Conservation Law
@add_constraints Demethy_TF_MC begin
  Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end
using Interact
@manipulate for KO = 0:0.01:1.0, K_nt = 0:0.01:1.0, Kd = 0:0.01:1.0, a1 = 0:0.1:10.0, d = 0:0.1:10.0, a_nt = 0:10:1000.0, aO = 0:10:1000.0, alphaT = 0:0.1:10.0, alphaO = 0:0.1:10.0, alphaN = 0:0.1:10.0, delta = 0:0.1:10.0, a_dn = 0:0.1:10.0, kh = 0:0.1:10.0, beta= 0:0.1:10.0
    p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta, 0., 0.05, 0.]
    ss = steady_states(Demethy_TF_MC,p)
    sort!(ss, by = x -> x[1])
    sb = stability(ss,Demethy_TF_MC,p)

    round_ss = [round.(i,digits = 2) for i in ss]
    dfc = DataFrame(vcat(round_ss))
    dfc.name = Demethy_TF_MC.syms
    var = [:N, :T, :O]
    df_iPS = dfc |> @filter(_.name in var) |> DataFrame
    @show Matrix(df_iPS)

    eigs = [round.(real(JE_stability(i, Demethy_TF_MC, p)[2]),digits = 2) for i in ss]
    @show eigs # show eigen values real part 2digits
    @show string.(sb)
    plot(sort([i[[1,2,9]] for i in ss]), xlabel = "N-T-O", label =string.(sb), marker = (:hexagon, 10, 0.7, stroke(1, 0.1, :black, :dot)))
end
