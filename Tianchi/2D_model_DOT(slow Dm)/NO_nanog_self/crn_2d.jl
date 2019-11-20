# ================ CRN ======================
using DiffEqBiological, LinearAlgebra, JLD2
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
using Interact;gr()
slider_rg = 0:10.0
@manipulate for KO=slider(0:0.01:1.0, value=0.3), K_nt = 0:0.01:1.0, Kd = 0:0.01:1.0, a1 = 0:0.1:10.0, d=0:0.1:10.0,  a_nt = 0:10:1000.0, aO = 0:10:1000.0, alphaT = 0:0.1:10.0, alphaO = 0:0.1:10.0, alphaN = 0:0.1:10.0, delta = 0:0.1:10.0, gamma = slider_rg, theta = 0:0.1:40.
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta, 0., 0.05, 0.] # version 1/2
    # p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta, 0., 0.05, 0.] # version 3
    ss = steady_states(Demethy_TF_MC,p)
    sort!(ss, by = x -> x[1])

    ss_round = [round.(i, digits = 3) for i in ss]
    dfc = DataFrame(vcat(ss_round))
    dfc.name = Demethy_TF_MC.syms
    # var = [:N, :T, :O, :Dm]
    var = Demethy_TF_MC.syms
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
    plot(sort([i[m2_idx] for i in ss]),ylims =(0,1),xticks = 1.:1.:3,label =string.(sb2),marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
    plot!(xticks = ([1.:1.:3;], ["N", "T", "O"]))
end






# We call this 2D model because we want to reduce the full model to only two genes with slow methylation dynamics (version 2 model).
# 2 Genes : oct4, Nanog
# The reduced model is actually 3D in terms of O,N, Dm
using ParameterizedFunctions
reduced_ODE_3d = @ode_def_bare begin # m1 -> N, m2 -> T,  m3 -> O, [N,O,Dm]
    # Asuume T also fast
    dN    = m1 - N*delta + (O*alphaN)/(KO + O)

    dO    = -(O^2*delta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - O*alphaO*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - KO*m3*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - O*m3*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*O*alphaO*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + KO*O*delta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) - N*O^2*alphaO*alphaT + N*O^3*alphaT*delta - N*O^2*alphaO*m2 - N*O^2*alphaT*m3 + N*O^3*delta*m2 - KO^2*N*m2*m3 - N*O^2*m2*m3 + K_nt*Kd*O^3*delta^2 - KO*N*O*alphaO*m2 - KO*N*O*alphaT*m3 - 2*KO*N*O*m2*m3 + Dm*N*O^2*alphaO*alphaT + K_nt*Kd*O^2*alphaO*delta + KO*N*O^2*alphaT*delta + Dm*N*O^2*alphaO*m2 - KO^2*K_nt*Kd*delta*m3 - K_nt*Kd*O^2*delta*m3 + 2*KO*N*O^2*delta*m2 + KO^2*N*O*delta*m2 + 2*KO*K_nt*Kd*O^2*delta^2 + KO^2*K_nt*Kd*O*delta^2 + KO*K_nt*Kd*O*alphaO*delta + Dm*KO*N*O*alphaO*m2 - 2*KO*K_nt*Kd*O*delta*m3 - Dm*K_nt*Kd*O^2*alphaO*delta - Dm*KO*K_nt*Kd*O*alphaO*delta)/((KO + O)*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + N*m2*(KO + O)^2 + N*O*alphaT*(KO + O) + K_nt*Kd*delta*(KO + O)^2)

    dDm   = -(Dm*N^2*O^2*alphaT^2*theta - 2*KO^2*K_nt*Kd^2*delta^2*gamma + Dm*KO^2*N^2*m2^2*theta + Dm*N^2*O^2*m2^2*theta - 2*KO*K_nt*Kd^2*O*delta^2*gamma + 2*Dm*KO*N^2*O*m2^2*theta + 2*Dm*N^2*O^2*alphaT*m2*theta + 2*Dm*KO^2*K_nt*Kd^2*delta^2*gamma + Dm*N*O*alphaT*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*KO*N*m2*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + Dm*N*O*m2*theta*(KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + 2*Dm*KO*N^2*O*alphaT*m2*theta + 2*Dm*KO*K_nt*Kd^2*O*delta^2*gamma - Dm*K_nt*Kd*N*O^2*alphaT*delta*theta + Dm*KO^2*K_nt*Kd*N*delta*m2*theta + Dm*K_nt*Kd*N*O^2*delta*m2*theta - Dm*KO*K_nt*Kd*N*O*alphaT*delta*theta + 2*Dm*KO*K_nt*Kd*N*O*delta*m2*theta)/(Kd*delta*(KO + O)*((KO^2*K_nt^2*Kd^2*delta^2 + 2*KO^2*K_nt*Kd*N*delta*m2 + KO^2*N^2*m2^2 + 2*KO*K_nt^2*Kd^2*O*delta^2 - 2*KO*K_nt*Kd*N*O*alphaT*delta + 4*KO*K_nt*Kd*N*O*delta*m2 + 2*KO*N^2*O*alphaT*m2 + 2*KO*N^2*O*m2^2 + K_nt^2*Kd^2*O^2*delta^2 - 2*K_nt*Kd*N*O^2*alphaT*delta + 2*K_nt*Kd*N*O^2*delta*m2 + N^2*O^2*alphaT^2 + 2*N^2*O^2*alphaT*m2 + N^2*O^2*m2^2)^(1/2) + N*O*alphaT + KO*N*m2 + N*O*m2 + KO*K_nt*Kd*delta + K_nt*Kd*O*delta))
end KO K_nt Kd a1 a_nt aO alphaT alphaO alphaN delta gamma theta m1 m2 m3


# 3 Genes: {N,T,O}, but with slow dynamics Dm. We could defind the BOA in this 4d dimensions
reduced_ODE_4d = @ode_def_bare begin
    dN  = m1 - N*delta + (O*alphaN)/(KO + O)
    dT  = m2 - T*delta + (N*O*T*alphaT)/((N*T + K_nt*Kd)*(KO + O))
    dO  = -(K_nt*Kd*O^2*delta + N*O^2*T*delta - KO*K_nt*Kd*m3 - N*O*T*alphaO - K_nt*Kd*O*m3 - KO*N*T*m3 - N*O*T*m3 + KO*N*O*T*delta + KO*K_nt*Kd*O*delta + Dm*N*O*T*alphaO)/((N*T + K_nt*Kd)*(KO + O))
    dDm = -(Dm*KO*K_nt*Kd^2*gamma - KO*K_nt*Kd^2*gamma + Dm*KO*N^2*T^2*theta + Dm*N^2*O*T^2*theta + Dm*KO*K_nt*Kd*N*T*theta + Dm*K_nt*Kd*N*O*T*theta)/(Kd*(N*T + K_nt*Kd)*(KO + O))
end KO K_nt Kd a1 a_nt aO alphaT alphaO alphaN delta gamma theta m1 m2 m3


# Test reduced 4d to determine tspan
# plotly()
# using DifferentialEquations
# u0 = [1.5, 0.0, 1.1666666666666667, 0.1111111111111111]
# tspan = (0., 5e2)
# p = [0.3, 0.2, 0.1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 0., 0.05, 0.]
# prob = ODEProblem(reduced_ODE_4d,u0,tspan,p)
# sol = solve(prob,Rosenbrock23())
# plot(sol)



# Add multithreading
using Base.Threads

# ====== Basin of Attraction (BOA) ========= testing now
# using Suppressor
# range = 0:10.
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
            NTO_Dm_idx = [1,2,8,15]
            Low = ss[1][NTO_Dm_idx]; Mid = ss[2][NTO_Dm_idx]; High = ss[3][NTO_Dm_idx]
            # @show Low Mid High
            # @show Low
            smpl_max = extrema(vcat([each[NTO_Dm_idx] for each in ss]...))[2]*1.5
            # @show smpl_max
            cube_O = cube_N = cube_T = LinRange(0.,smpl_max,10)
            con = 1
            tspan = (0., 1e2)
            soma = []
            TV = 0
            for O = cube_O, N = cube_N, T = cube_T, Dm = LinRange(0., con, 6)
                u0 = [N,T,O,Dm]
                # @show u0
                prob = ODEProblem(reduced_ODE_4d,u0,tspan,p)
                sol = solve(prob,SSRootfind())
                f_ss = norm(sol[end] .- Low) < 0.1 ? 1 : 0
                # @show f_ss
                push!(soma, f_ss)
                TV += 1
            end
            DOT = sum(soma)/TV
            C[i] = DOT
        else
            C[i] = NaN
        end
    end
    return C
end

p = [0.3, 0.2, 0.1, 1., 100., 100., .8, 1.0, 1.0, 1., 5., 5.] # version 2
# rg = [0:10.; 10.:10.:50.]
rg = 0:0.8:10.
@time C = DOT_Volume_params(p..., range = rg,param = "γ")
plot(rg, C, ylabel = "Basin of Attraction",xlabel = "Gamma", legend = false)

# get non-NAN values index and plot in log10 scale
v_idx = @. !isnan.(C)
plot(rg[v_idx], C[v_idx], xscale = :log10,ylabel = "Basin of Attraction",xlabel = "Theta", legend = false)







# BOA plots
using Distributed;addprocs(4)
using ProgressMeter
# pv = plot();
Cγ = []; Cθ = []; rg = exp10.(-2:0.1:2.5) #0:5.
# rg = [0:10.; 10.:10.:200.]
@showprogress "Computing..." for KO = 0.3, K_nt = 0.2, Kd = 0.1, a1 = 1., a_nt = 100., aO = 100., alphaT = 1., alphaO = 1., alphaN = 1., delta = 1:5, gamma = 1:2.:20., theta = 1.
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta]
    C = DOT_Volume_params(p..., range = rg,param = "θ")
    # @show C
    # push!(Cγ,C) #  Methylation
    push!(Cθ,C)  # Demethylation
    sleep(0.01)
end

plot(rg, Cθ.*100, ylabel = "Basin of Attraction %", xscale =:log10, xlabel = "Theta" ,legend =false)




# ======= γ pplot ======================
Cγ =[];rg = exp10.(-2:0.1:2.5)
D_rg = 0.1:0.1:0.5
@time @showprogress "Computing..." for KO = D_rg, K_nt = D_rg, Kd = D_rg, a1 = 1., a_nt = 100., aO = 100., alphaT = 1., alphaO = 1., alphaN = 1., delta = 1:5, gamma = 10., theta = 1:2.:20.
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta]
    C = DOT_Volume_params(p..., range = rg, param = "γ")
    # push!(Cm,C)
    push!(Cγ,C)
    sleep(0.01)
end

# below just test for a potential lower bound
@time @showprogress "Computing..." for KO = 0.3, K_nt = 0.2, Kd = 0.1, a1 = 10., a_nt = 110., aO = 500., alphaT = 5., alphaO = 5., alphaN = 5., delta = 1, gamma = 10., theta = 5.
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta]
    C = DOT_Volume_params(p..., range = rg, param = "γ")
    # push!(Cm,C)
    push!(Cγ,C)
    sleep(0.01)
end
# @save "BOA_γ[KO, K_nt,Kd, delta, theta]{rg = exp10.(-2:0.1:2.5),Cγ }.jld2" Cγ
@load "BOA_γ[KO, K_nt,Kd, delta, theta]{rg = exp10.(-2:0.1:2.5),Cγ }.jld2" Cγ
# BOA_γ_plt = plot(rg, Cγ, ylabel = "Basin of Attraction", xlabel = "Gamma" ,legend =false) # regulat scale
# BOA_γ_plt = plot(rg, Cγ, ylabel = "Basin of Attraction", xscale =:log10, xlabel = "gamma" ,legend =false) # log10 scale
BOA_γ_plt = plot(rg, Cγ.*100, ylabel = "Basin of Attraction %", xscale =:log10, xlabel = L"\gamma" ,legend =false) # ylabel
savefig(BOA_γ_plt, "~/Desktop/BOAp_γ_full_model.png")



# ======= θ plots =================================
Cθ =[];rg = exp10.(-2:0.1:2.5)
D_rg = 0.1:0.1:0.5
@time @showprogress "Computing..." for KO = D_rg, K_nt = D_rg, Kd = D_rg, a1 = 1., a_nt = 100., aO = 100., alphaT = 1., alphaO = 1., alphaN = 1., delta = 1:5, gamma = 1:2.:20, theta = 1.
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta]
    C = DOT_Volume_params(p..., range = rg, param = "θ")
    # push!(Cm,C)
    push!(Cθ,C)
    sleep(0.01)
end
@save "BOA_θ[KO,K_nt,Kd,delta,γ]{rg = exp10.(-2:0.1:2.5),Cθ}.jld2" Cθ
@load "BOA_θ[KO,K_nt,Kd,delta,γ]{rg = exp10.(-2:0.1:2.5),Cθ}.jld2" Cθ
BOA_θ_plt = plot(rg, Cθ.*100, ylabel = "Basin of Attraction %", xscale =:log10, xlabel = L"\theta" ,legend =false)
savefig(BOA_θ_plt, "~/Desktop/BOAp_θ_full_model.png")



# ======== Keep R = γ/θ same and change γ =================
function BOA_Volume_params_R_fixed(KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, R; range = 0:0.1:10., param = "γ")
    C = similar(range)
    for i in eachindex(range) #@suppress
        # @show i
        if param == "γ"
            # Still changing γ,but fixing R
            p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, range[i], range[i]/R, 0., 0.05, 0.]
        end
        ss = steady_states(Demethy_TF_MC,p)
        sort!(ss, by = x -> x[1])
        # @show ss[1]
        if length(ss) >2
            NTO_Dm_idx = [1,2,8,15]
            Low = ss[1][NTO_Dm_idx]; Mid = ss[2][NTO_Dm_idx]; High = ss[3][NTO_Dm_idx]
            # @show Low Mid High
            # @show Low
            smpl_max = extrema(vcat([each[NTO_Dm_idx] for each in ss]...))[2]*1.5
            # @show smpl_max
            cube_O = cube_N = cube_T = LinRange(0.,smpl_max,10)
            con = 1
            tspan = (0., 1e2)
            soma = []
            TV = 0
            for O = cube_O, N = cube_N, T = cube_T, Dm = LinRange(0., con, 6)
                u0 = [N,T,O,Dm]
                # @show u0
                prob = ODEProblem(reduced_ODE_4d,u0,tspan,p)
                sol = solve(prob,SSRootfind())
                f_ss = norm(sol[end] .- Low) < 0.1 ? 1 : 0
                # @show f_ss
                push!(soma, f_ss)
                TV += 1
            end
            DOT = sum(soma)/TV
            C[i] = DOT
        else
            C[i] = NaN
        end
    end
    return C
end

CRγ =[];rg = exp10.([-3:0.2:0; 0:0.2:2.5])
# D_rg = 0.1:0.1:0.5
@time @showprogress "Computing..." for KO = 0.3, K_nt = 0.2, Kd = 0.1, a1 = 1., a_nt = 100., aO = 100., alphaT = 1., alphaO = 1., alphaN = 1., delta = 1., R = [0.01, 0.1, 1., 10.,20.,30.], gamma = 1.
    pγ = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, R]
    Cγ = BOA_Volume_params_R_fixed(pγ..., range = rg, param = "γ")
    push!(CRγ,Cγ)
    # pθ = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, theta*R, theta]
    # Cθ = DOT_Volume_params(pθ..., range = rg, param = "θ")
    # push!(CRθ,Cθ)
    sleep(0.01)
end
@save "BOA_R_γ[R 5orders of mag]{rg = exp10.([-3:0.2:0; 0:0.2:2.5]), CRγ}.jld2" CRγ
@load "BOA_R_γ[R 5orders of mag]{rg = exp10.([-3:0.2:0; 0:0.2:2.5]), CRγ}.jld2" CRγ
BOA_R_γ_plt = plot(rg, CRγ.*100, ylabel = "Basin of Attraction %", xscale =:log10, xlabel = L"\gamma" , label =["0.01", "0.1" , "1" , "10", "20", "30"], legendtitle = "R ratio", legend = :topright)
BOA_R_γ_plt = plot(rg[5:end], [i[5:end].*100 for i in CRγ], ylabel = "Basin of Attraction %", xscale =:log10, xlabel = L"\gamma" , label =["0.01", "0.1" , "1" , "10", "20", "30"], legendtitle = "R ratio", legend = :topright) # ignore the samll rg region
# BOA_R_θ_plt = plot(rg, CRθ.*100, ylabel = "Basin of Attraction %", xscale =:log10, xlabel = "Theta" , label =["0.01", "0.1" , "1" , "10", "20", "30"], legendtitle = "R ratio", legend = :topright)
# BOA_R_plt = plot(BOA_R_γ_plt,BOA_R_θ_plt, layout = (1,2))
savefig(BOA_R_γ_plt, "~/Desktop/BOAp_R_4orders.png")



@time @showprogress "Computing..." for KO = 0.5, K_nt = 0.2, Kd = 0.1, a1 = 1., a_nt = 100., aO = 100., alphaT = 1., alphaO = 1., alphaN = 1., delta = 1., R = [0.01, 0.1, 1., 10.,20.,30.], gamma = 1.
    pγ = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, R]
    Cγ = BOA_Volume_params_R_fixed(pγ..., range = rg, param = "γ")
    push!(CRγ,Cγ)
    # pθ = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, theta*R, theta]
    # Cθ = DOT_Volume_params(pθ..., range = rg, param = "θ")
    # push!(CRθ,Cθ)
    sleep(0.01)
end

BOA_R_γ_plt_KO = plot(rg, CRγ.*100, ylabel = "Basin of Attraction %", xscale =:log10, xlabel = "Gamma" , label =["0.01", "0.1" , "1" , "10", "20", "30"], legendtitle = "R ratio", legend = :topright)
