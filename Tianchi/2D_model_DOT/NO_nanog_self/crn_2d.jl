# ================ CRN ======================
using DiffEqBiological, LinearAlgebra
using Plots;gr()
using DataFrames, Queryverse, Latexify
include(pwd()*"/functions.jl")

# --- CRN model ------

# === version 1.
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT  #💚 NO dimerization

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

    # --- Protein production
    alphaT,                   Dt11 → Dt11 + T
    alphaO,                   Do11 → Do11 + O
    alphaN,                   Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),           (N, T, O) → ∅

    # ============   🔺 ===============
    # # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    # a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    # kh,                  NT + D5mc → D5hmc + NT     # NT oxidize 5mc -> 5hmC
    # beta,                D5hmc →  Do00              # 5hmC -> C by
    # ============ ⟺ ===================
    (gamma, theta),           Do00 ↔ Dm # simplification of the 🔺

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end KO K_nt Kd a1 a_nt aO alphaT alphaO alphaN delta gamma theta m1 m2 m3

# === version 2. NT control the demethylation
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
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

    # --- Protein production
    alphaT,                   Dt11 → Dt11 + T
    alphaO,                   Do11 → Do11 + O
    alphaN,                   Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (delta,delta,delta),           (N, T, O) → ∅

    # ============   🔺 ===============
    # # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    # a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    # kh,                  NT + D5mc → D5hmc + NT     # NT oxidize 5mc -> 5hmC
    # beta,                D5hmc →  Do00              # 5hmC -> C by
    # ============ ⟺ ===================
    gamma,             Do00 → Dm   # simplification of the 🔺
    theta,             Dm + NT ⇀ Do00 + NT
    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end KO K_nt Kd a1 a_nt aO alphaT alphaO alphaN delta gamma theta m1 m2 m3

# === version 3. NT2 control the demethylation and promoters [prefered]
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    (d,d),                    NT + NT ↔ NT2 # NT dimerization
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
    # # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    # a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    # kh,                  NT + D5mc → D5hmc + NT     # NT oxidize 5mc -> 5hmC
    # beta,                D5hmc →  Do00              # 5hmC -> C by
    # ============ ⟺ ===================
    gamma,             Do00 → Dm   # simplification of the 🔺
    theta,             Dm + NT2 ⇀ Do00 + NT2
    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta gamma theta m1 m2 m3

# Add constraints
@add_constraints Demethy_TF_MC begin
  # Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1 # if use 🔺for demethylation
  Do00 + Do01 + Do10 + Do11 + Dm  = 1 # if use ⟺ for dedemethylation
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end




p = [0.3, 0.2, 0.1, 1., 10., 500., 500., 1.0, 1.0, 1.0, 1., 1., 1.,  0., 0.05, 0.] # version 3
p = [0.3, 0.2, 0.1, 1., 500., 500., 1.0, 1.0, 1.0, 1., 1., 1.,  0., 0.05, 0.] # version 2
ss = steady_states(Demethy_TF_MC,p)
sort!(ss, by = x -> x[1])

sb = stability(ss,Demethy_TF_MC,p)
sb2 = stability_tianchi(ss,Demethy_TF_MC,p,3)



dfc = DataFrame(vcat(ss))
dfc.name = Demethy_TF_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame




# First Visulization =====================
# ===== 3d model DOT defined by N-T-O-----
using Interact
slider_rg = 0:10.0
@manipulate for KO = 0:0.01:1.0, K_nt = 0:0.01:1.0, Kd = 0:0.01:1.0, a1 = 0:0.1:10.0, d=0:0.1:10.0,  a_nt = 0:10:1000.0, aO = 0:10:1000.0, alphaT = 0:0.1:10.0, alphaO = 0:0.1:10.0, alphaN = 0:0.1:10.0, delta = 0:0.1:10.0, gamma = slider_rg, theta = slider_rg
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta, 0., 0.05, 0.] # version 1/2
    # p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta, 0., 0.05, 0.] # version 3
    ss = steady_states(Demethy_TF_MC,p)
    sort!(ss, by = x -> x[1])

    # ss_round = [round.(i, digits = 3) for i in ss]
    dfc = DataFrame(vcat(ss))
    dfc.name = Demethy_TF_MC.syms
    var = [:N, :T, :O]
    @show dfc1 = dfc |> @filter(_.name in var) |> DataFrame

    if length(ss) >2
        DOT   = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3])) # 2nd def
        @show DOT
    else
         @show "single state"
    end

    # Check stability
    sb2 = stability_tianchi(ss,Demethy_TF_MC,p,3)

    # Plotting
    m2_idx = [1,2,8]
    plot(sort([i[m2_idx] for i in ss]),xticks = 1.:1.:3,label =string.(sb2),marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
    plot!(xticks = ([1.:1.:3;], ["N", "T", "O"]))
end





# We call this 2D model because we want to reduce the full model to only two genes with slow methylation dynamics (version 2 model).
# 2 Genes : oct4, Nanog
# The reduced model is actually 4D in terms of O,N, Do00, Dm
using ParameterizedFunctions
reduced_ODE_4d = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN    = m1 - N*delta + (O*alphaN)/(KO + O)

    dDo00 = (Dm*N^2*O^2*alphaT^2*theta - 2*KO^2*K_nt*Kd^2*delta^2*gamma + Dm*KO^2*N^2*m2^2*theta + Dm*N^2*O^2*m2^2*theta - 2*KO*K_nt*Kd^2*O*delta^2*gamma + 2*Dm*KO*N^2*O*m2^2*theta + 2*Dm*N^2*O^2*alphaT*m2*theta + 2*Dm*KO^2*K_nt*Kd^2*delta^2*gamma + Dm*N*O*alphaT*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*KO*N*m2*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*N*O*m2*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + 2*Dm*KO*N^2*O*alphaT*m2*theta + 2*Dm*KO*K_nt*Kd^2*O*delta^2*gamma - Dm*K_nt*Kd*N*O^2*alphaT*delta*theta + Dm*KO^2*K_nt*Kd*N*delta*m2*theta + Dm*K_nt*Kd*N*O^2*delta*m2*theta - Dm*KO*K_nt*Kd*N*O*alphaT*delta*theta + 2*Dm*KO*K_nt*Kd*N*O*delta*m2*theta)/(Kd*delta*(KO + O)*((KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + N*O*alphaT + KO*N*m2 + N*O*m2 + KO*K_nt*Kd*delta + K_nt*Kd*O*delta))

    dO    = -(O^2*delta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - O*alphaO*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - KO*m3*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - O*m3*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*O*alphaO*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + KO*O*delta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - N*O^2*alphaO*alphaT + N*O^3*alphaT*delta - N*O^2*alphaO*m2 - N*O^2*alphaT*m3 + N*O^3*delta*m2 - KO^2*N*m2*m3 - N*O^2*m2*m3 + K_nt*Kd*O^3*delta^2 - KO*N*O*alphaO*m2 - KO*N*O*alphaT*m3 - 2*KO*N*O*m2*m3 + Dm*N*O^2*alphaO*alphaT + K_nt*Kd*O^2*alphaO*delta + KO*N*O^2*alphaT*delta + Dm*N*O^2*alphaO*m2 - KO^2*K_nt*Kd*delta*m3 - K_nt*Kd*O^2*delta*m3 + 2*KO*N*O^2*delta*m2 + KO^2*N*O*delta*m2 + 2*KO*K_nt*Kd*O^2*delta^2 + KO^2*K_nt*Kd*O*delta^2 + KO*K_nt*Kd*O*alphaO*delta + Dm*KO*N*O*alphaO*m2 - 2*KO*K_nt*Kd*O*delta*m3 - Dm*K_nt*Kd*O^2*alphaO*delta - Dm*KO*K_nt*Kd*O*alphaO*delta)/((KO + O)*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + N*m2*(KO + O)^2 + N*O*alphaT*(KO + O) + K_nt*Kd*delta*(KO + O)^2)

    dDm   = -(Dm*N^2*O^2*alphaT^2*theta - 2*KO^2*K_nt*Kd^2*delta^2*gamma + Dm*KO^2*N^2*m2^2*theta + Dm*N^2*O^2*m2^2*theta - 2*KO*K_nt*Kd^2*O*delta^2*gamma + 2*Dm*KO*N^2*O*m2^2*theta + 2*Dm*N^2*O^2*alphaT*m2*theta + 2*Dm*KO^2*K_nt*Kd^2*delta^2*gamma + Dm*N*O*alphaT*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*KO*N*m2*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*N*O*m2*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + 2*Dm*KO*N^2*O*alphaT*m2*theta + 2*Dm*KO*K_nt*Kd^2*O*delta^2*gamma - Dm*K_nt*Kd*N*O^2*alphaT*delta*theta + Dm*KO^2*K_nt*Kd*N*delta*m2*theta + Dm*K_nt*Kd*N*O^2*delta*m2*theta - Dm*KO*K_nt*Kd*N*O*alphaT*delta*theta + 2*Dm*KO*K_nt*Kd*N*O*delta*m2*theta)/(Kd*delta*(KO + O)*((KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + N*O*alphaT + KO*N*m2 + N*O*m2 + KO*K_nt*Kd*delta + K_nt*Kd*O*delta))
end KO K_nt Kd a1 a_nt aO alphaT alphaO alphaN delta gamma theta m1 m2 m3



# ====== Basin of Attraction (BOA) ========= testing now
using Suppressor
range = 0:10.
function DOT_Volume_params(KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta; range = 0:0.1:10., param = "γ")
    C = similar(range)
    for i in eachindex(range) #@suppress
        # @show i
        if param == "γ"
            p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, range[i], theta, 0., 0.05, 0.]
        elseif param == "θ"
            p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, range[i],  0., 0.05, 0.]
        end
        ss = steady_states(Demethy_TF_MC,p)
        sort!(ss, by = x -> x[1])
        # @show ss[1]
        if length(ss) >2
            N_Do_O_Dm_idx = [1,4,8,15]
            Low = ss[1][N_Do_O_Dm_idx]; Mid = ss[2][N_Do_O_Dm_idx]; High = ss[3][N_Do_O_Dm_idx]
            @show Low Mid High

            smpl_max = extrema(vcat([each[N_Do_O_Dm_idx] for each in ss]...))[2]*1.5
            @show smpl_max
            cube_O = cube_N = LinRange(0.,smpl_max,15)
            con = 1
            tspan = (0., 5e2)
            soma = []
            TV = 0
            for O = cube_O, N = cube_N, Do00 = LinRange(0., con, 10)
                if con - Do00 >0
                    Dm = con - Do00
                    u0 = [N,Do00,O,Dm]
                    # @show u0
                    prob = ODEProblem(reduced_ODE_4d,u0,tspan,p)
                    sol = solve(prob,Rosenbrock23())
                    rd_idx = [1,4,8,15]
                    f_ss = norm(sol[end] .- ss[1][rd_idx]) < 0.1 ? 1 : 0
                    push!(soma, f_ss)
                    TV += 1
                    # @show f_ss
                end
            end
            DOT = sum(soma)/TV
            C[i] = DOT
        else
            C[i] = NaN
        end
    end
    return C
end

p = [0.3, 0.2, 0.1, 1., 10., 10., .5, 1.0, 1.0, 1., 10., 1.] # version 2
rg = [0:10.; 10.:10.:50.]
C = DOT_Volume_params(p..., range = rg,param = "θ")
plot(rg, C)




# Test reduced 4d
using DifferentialEquations
u0 = [5.698513885281152, 0.2631578947368421, 3.256293648732086, 4.7368421052631575]
tspan = (0., 5e2)
p = [0.3, 0.2, 0.1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 0., 0.05, 0.]
prob = ODEProblem(reduced_ODE_4d,u0,tspan,p)
sol = solve(prob,Rosenbrock23())
plot(sol)
