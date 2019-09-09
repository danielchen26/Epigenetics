using DiffEqBiological, DifferentialEquations
using DataFrames, LinearAlgebra, Plots
using Queryverse

module epi
    include("../functions.jl")
end



# --------Paper model -------
ps = @reaction_network begin
    # N T complex
    (a₁, Kd*a₁),         N + T ↔ NT
    (d, d),             NT + NT ↔ NT2

    # --- Promoter binding
    # Oct4
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₀ᴼ ↔ D₁₀ᴼ
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₁ᴼ ↔ D₁₁ᴼ
    (aₒ, Kₒ*aₒ),          O + D₀₀ᴼ ↔ D₀₁ᴼ
    (aₒ, Kₒ*aₒ),          O + D₁₀ᴼ ↔ D₁₁ᴼ
    # TET
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₀ᵀ ↔ D₁₀ᵀ
    (aₙₜ, Kₙₜ*aₙₜ),        NT2 + D₀₁ᵀ ↔ D₁₁ᵀ
    (aₒ, Kₒ*aₒ),          O + D₀₀ᵀ ↔ D₀₁ᵀ
    (aₒ, Kₒ*aₒ),          O + D₁₀ᵀ ↔ D₁₁ᵀ
    # Nanog
    (aₒ, Kₒ*aₒ),          O + D₀ᴺ ↔ D₁ᴺ

    # --- Protein production
    αₜ,                 D₁₁ᵀ → D₁₁ᵀ + T
    αₒ,                 D₁₁ᴼ → D₁₁ᴼ + O
    αₙ,                 D₁ᴺ → D₁ᴺ + N

    # --- Dilution and Degradation
    (δ,δ,δ),           (N, T, O) → ∅
    #(β,β),             (NT, NT2) → ∅

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ m1 m2 m3

p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1,   0., 0.05, 0.]
@add_constraints ps begin
  D₀₀ᴼ + D₁₀ᴼ + D₀₁ᴼ + D₁₁ᴼ = 1
  D₀₀ᵀ + D₁₀ᵀ + D₀₁ᵀ + D₁₁ᵀ = 1
  D₀ᴺ  + D₁ᴺ = 1
end
ss2 = steady_states(ps,p)
stability(ss2,ps,p)

df = DataFrame(vcat(ss2))
df.name = ps.syms
var = [:N, :T, :O]
df2 = df |> @filter(_.name in var) |> DataFrame

# N_sample = 2
# for u01 in epi.randfixsum(N_sample,4,1) , u02 in epi.randfixsum(N_sample,4,1), u03 in epi.randfixsum(N_sample,2,1)
#     u0 = [rand(1,4) u01 rand(1) u02 u03]
#     tspan = (0,1e6)
#     oprob = ODEProblem(ps, u0, tspan, p)
#     osol  = solve(oprob, Rosenbrock23())
#     # plot(osol,  ylims = (0,2))
#
#     # calculate eig_spectrum
#     t = 0.0
#     out = rand(15,15)
#     ps.jac(out,osol[end],p,t)
#     # @show out
#     eig = eigen(out).values
#     println("Eigen value:\n")
#     display(eig)
# end
# gif(anim, fps = 1)





# =========############# =========#############===============#############==========
# =========############# =========#############===============#############==========




#   ------ Methylation CRN ------Matlab_syntax:NC
Demethy_crn = @reaction_network begin
    # N T complex
    (a₁, Kd*a₁),       N + T ↔ NT        # T is recruited by N

    # --- Promoter binding
    # Oct4 cycle
    (r₁, Kₒ*r₁ ),      O + Dₒ ↔ Dᵒ           # O Auto-activation
    α₁,               Dᵒ → Dᵒ + O           # O translation
    δₒ,               O → ∅                 # Deg
    # Oct4 de-Methylation cycle
    a0,               Dₒ → D̄               # D + DNMT ↔ C₁ → D̄ + DNMT
    β,                Dₕ →  Dₒ               # 5hmC -> C by dilution(replication)
    # TET protein
    (r₁,Kₒ*r₁),       Dₜ + O ↔ Dᵗ            # O -> T promoter
    α₁,               Dᵗ → Dᵗ + T            # T translation
    δₒ,               T → ∅                 # Deg
    # Nanog cycle
    (r₁, Kₒ*r₁ ),       N + Dₙ ↔ Dᴺ            # N Auto-activation
    α₁,                Dᴺ → Dᴺ + N            # N translation
    δₒ,                N → ∅                  # Deg
    # O activate N
    (ζ₁, Kₒ*ζ₁),       O + Dₙ ↔ Dᴺ            # O -> N promoter

    # ---- Nanog guided Tet1/2
    k₃,               NT + D̄ → Dₕ + NT       # NT oxidize 5mc -> 5hmC
    # k₄,               NT + Dₜ → Dᵗ + NT       # NT induce not only Oct4 but T1

    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O

end Kₒ Kd r₁ a₁ α₁ δₒ a0 β ζ₁ k₃ m1 m2 m3# k₄

@add_constraints Demethy_crn begin
  Dₒ + Dᵒ + D̄ + Dₕ  = 1
  Dₜ + Dᵗ = 1
  Dₙ + Dᴺ = 1
end

params = [0.3, 0.1, 1, 1, 1, 1, 0.5, 1, 1, 1,   0., 0.05, 0.]
ss = steady_states(Demethy_crn,params)
stability(ss,Demethy_crn,params)

df1 = DataFrame(vcat(ss))
df1.name = Demethy_crn.syms
var = [:N, :T, :O]
@show df11 = df1 |> @filter(_.name in var) |> DataFrame

bif_grid_dia = bifurcation_grid_diagram(Demethy_crn, params, :ζ₁, 1.:10., :k₃, (1, 5))
plot(bif_grid_dia)

plotly()





#  explortion of params
for i  = exp10.(1:.1:8)
    params[14:19] .= i
    # @show ss = steady_states(O_N_auto,params)
    dfi = DataFrame(vcat(ss))
    dfi.name = O_N_auto.syms
    var = [:N, :T, :O]
    dff = dfi |> @filter(_.name in var) |> DataFrame
    df_final = unstack(stack(dff),:variable,:name,:value)
    df_final.stability = stability(ss,O_N_auto,params)
    @show df_final
end






# =========############# =========#############===============#############==========
# =========############# =========#############===============#############==========

using DiffEqBiological
using DataFrames, Queryverse
using Plots; plotly()

#   ------ Methylation CRN ------Matlab_syntax:C

Demethy_crn_MatlabC = @reaction_network begin
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


@add_constraints Demethy_crn_MatlabC begin
  Do + DO + D5mc + D5hmc  = 1
  Dt + DT = 1
  Dn + DN = 1
end

params = [0.3, 0.1, 1.63792, 1.20623, 1.28428, 1.28329, 0.0776827, 1.46203, 1.75094, 1.51652,    0., 0.05, 0. ]
ss = steady_states(Demethy_crn_MatlabC,params)
stability(ss,Demethy_crn_MatlabC,params)

function my_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,p,t)
    return (jac,eigen(jac).values)
end


J = rand(length(Demethy_crn_MatlabC.syms),length(Demethy_crn_MatlabC.syms) );t=0
jacfun(Demethy_crn_MatlabC)(J, ss[end], params, t)
eigen(J).values


J = rand(length(rn.syms),length(rn.syms) );t=0
jacfun(rn)(J, ss[end], params, t)
eigen(J).values



df1 = DataFrame(vcat(ss))
df1.name = Demethy_crn_MatlabC.syms
var = [:N, :T, :O]
@show df11 = df1 |> @filter(_.name in var) |> DataFrame

bif_grid_dia = bifurcation_grid_diagram(Demethy_crn_MatlabC, params, :K0, 0.:0.1:1.,  :Kd, (0.1,10.))

plot(bif_grid_dia)

#  explortion of params
for i = 1:20, j =  1:20
    params[10] = i
    params[9] = j
    ss = steady_states(Demethy_crn_MatlabC,params)
    dfi = DataFrame(vcat(ss))
    dfi.name = Demethy_crn_MatlabC.syms
    var = [:N, :T, :O]
    dff = dfi |> @filter(_.name in var) |> DataFrame
    df_final = unstack(stack(dff),:variable,:name,:value)
    df_final.stability = stability(ss,Demethy_crn_MatlabC,params)
    @show df_final
end






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


params = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
ss = steady_states(Demethy_TF_MC,params)
stability(ss,Demethy_TF_MC,params)

dfc = DataFrame(vcat(ss))
dfc.name = Demethy_TF_MC.syms
var = [:N, :T, :O]
@show dfc1 = dfc |> @filter(_.name in var) |> DataFrame
















# -------------------- doc for bifurcation_grid_diagram ---------
# color:
#Blue: All eigenvalues smaller than 0, and no imaginary parts.
# Cyan: All eigenvalues smaller than 0, and at least one with imaginary parts.
# Orange: At least one eigenvalues larger than 0, and at least one with imaginary parts.
# Red: At least one eigenvalues larger than 0, and no imaginary parts.
