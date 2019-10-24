# ================ CRN ======================
using DiffEqBiological, LinearAlgebra
using Plots;plotly()
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


# p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
p = [0.1, 0.2, 0.2, 5., 5., 500., 500., 5.0, 5.0, 5.0, 5., 5., 5., 5.,  0., 0.05, 0.]
@add_constraints Demethy_TF_MC begin
  Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end

ss = steady_states(Demethy_TF_MC,p)
stability(ss,Demethy_TF_MC,p)

function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,p,t)
    return (jac,eigen(jac).values)
end

@show JE_stability(ss[3], Demethy_TF_MC, p)[2]





















# Bifurcation diagram
bif_grid_dia = bifurcation_grid_diagram(Demethy_TF_MC, p, :KO, 0.:0.1:0.55,  :Kd, (0.1,.5))
bp = plot(bif_grid_dia)
plot!(bp,[[],[],[],[]],color=[:blue :cyan :orange :red],label = ["Stable Real" "Stable Complex" "Unstable Complex" "Unstable Real"])




dfc = DataFrame(vcat(ss))
dfc.name = Demethy_TF_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame













# First Visulization =====================
# ===== 3d model DOT defined by N-T-O-----

using Interact, Suppressor, ProgressMeter
@manipulate for KO = 0:0.01:1.0, K_nt = 0:0.01:1.0, Kd = 0:0.01:1.0, a1 = 0:0.1:10.0, d = 0:0.1:10.0, a_nt = 0:10:1000.0, aO = 0:10:1000.0, alphaT = 0:0.1:10.0, alphaO = 0:0.1:10.0, alphaN = 0:0.1:10.0, delta = 0:0.1:10.0, a_dn = 0:0.1:10.0, kh = 0:0.1:10.0, beta= 0:0.1:10.0
    p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta, 0., 0.05, 0.]
    ss = steady_states(Demethy_TF_MC,p)
    sort!(ss, by = x -> x[1])

    dfc = DataFrame(vcat(ss))
    dfc.name = Demethy_TF_MC.syms
    var = [:N, :T, :O]
    df_iPS = dfc |> @filter(_.name in var) |> DataFrame
    @show Matrix(df_iPS)
    if length(ss) >2
        DOT   = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3]))
        @show DOT
    else
         @show "single state"
    end
    plot(sort([i[[1,2,9]] for i in ss]), xlabel = "N-T-O",marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
    # println("The SS states are:")
    # sort!(ss, by = x -> x[1])
    # @show ss
    # @show [round.(EigenVector(rn1d,p,i),digits = 2) for i in ss]

    # eigs = [JE_stability(i,Demethy_TF_MC,p)[2] for i in ss]
    # println("The eigs are:")
    # # @show eigs
    # println("Eig values:",hcat(eigs...))
end
