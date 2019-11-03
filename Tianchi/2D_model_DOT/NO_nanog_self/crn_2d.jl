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

# === version 2. NT control the demethylation [prefered]
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
@add_constraints Demethy_TF_MC begin
  # Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1 # if use 🔺for demethylation
  Do00 + Do01 + Do10 + Do11 + Dm  = 1 # if use ⟺ for dedemethylation
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end

p = [0.3, 0.2, 0.1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1,  0., 0.05, 0.]
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
@manipulate for KO = 0:0.01:1.0, K_nt = 0:0.01:1.0, Kd = 0:0.01:1.0, a1 = 0:0.1:10.0,  a_nt = 0:10:1000.0, aO = 0:10:1000.0, alphaT = 0:0.1:10.0, alphaO = 0:0.1:10.0, alphaN = 0:0.1:10.0, delta = 0:0.1:10.0, gamma = 0:0.1:10.0, theta = 0:0.1:10.0
    p = [KO, K_nt, Kd, a1, a_nt, aO, alphaT, alphaO, alphaN, delta, gamma, theta, 0., 0.05, 0.]
    ss = steady_states(Demethy_TF_MC,p)
    sort!(ss, by = x -> x[1])

    ss_round = [round.(i, digits = 3) for i in ss]
    dfc = DataFrame(vcat(ss_round))
    dfc.name = Demethy_TF_MC.syms
    var = [:N, :T, :O]
    df_iPS = dfc |> @filter(_.name in var) |> DataFrame
    @show Matrix(df_iPS)

    if length(ss) >2
        DOT   = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3])) # 2nd def
        @show DOT
    else
         @show "single state"
    end

    # Check stability
    sb2 = stability_tianchi(ss,Demethy_TF_MC,p,3)

    # Plotting
    plot(sort([i[[1,2,9]] for i in ss]),xticks = 1.:1.:3,label =string.(sb2),marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
    plot!(xticks = ([1.:1.:3;], ["N", "T", "O"]))
end





# ====== Basin of Attraction (BOA) ========= Not yet started
range = 0:0.1:10.
function DOT_Volume_params(α, ϵ, β, η, γ, θ, δ; range = 0:0.1:10., param = "γ")
    C = similar(range)
    @suppress for i in eachindex(range)
        # @show i
        if param == "γ"
            p = [α, ϵ, β, η, range[i], θ, δ]
        elseif param == "θ"
            p = [α, ϵ, β, η, γ, range[i], δ]
        end
        ss = steady_states(rn1d_leak,p)
        sort!(ss, by = x -> x[1])
        # @show ss[1]
        if length(ss) >2
            smpl_max = extrema(vcat(ss...))[2]*1.5
            cube_O = LinRange(0.,smpl_max,10)
            con = 5
            tspan = (0., 5e2)
            soma = []
            TV = 0
            for O = cube_O, Di = LinRange(0., con, 20), Da = LinRange(0.,con,20)
                if con - Di -Da >0
                    Dm = con - Di -Da
                    u0 = [O,Di,Da,Dm]
                    # @show u0
                    prob = ODEProblem(rn1d_leak,u0,tspan,p)
                    sol = solve(prob,Rosenbrock23())
                    f_ss = norm(sol[end] .- ss[1]) < 0.1 ? 1 : 0
                    push!(soma, f_ss)
                    TV += 1
                    # @show sol[end]
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
# ===test above two functions ====
C_1d, C_3d = DOT_γ(5., 5., 5., 0.1, 5., 5.)
C = DOT_Volume_params(5., 5., 5., 0.1, 5., 0., 5., range = [0:0.1:2.;2.:1.:10.], param = "θ")
plot([0:0.1:2.;2.:1.:10.],C)






# dN     = m1 + Dn1*alphaN - N*delta - N*T*a1 + NT*Kd*a1
# dT     = m2 + Dt11*alphaT - T*delta - N*T*a1 + NT*Kd*a1
# dNT    = K_nt*a_nt*Do10 + K_nt*a_nt*Do11 + K_nt*a_nt*Dt10 + K_nt*a_nt*Dt11 + N*T*a1 - NT*Kd*a1 - NT*a_nt*Do00 - NT*a_nt*Do01 - NT*a_nt*Dt00 - NT*a_nt*Dt01
# dDo00  = -Do00*gamma + KO*aO*Do01 + K_nt*a_nt*Do10 + NT*Dm*theta - NT*a_nt*Do00 - O*aO*Do00
# dDo10  = KO*aO*Do11 - K_nt*a_nt*Do10 + NT*a_nt*Do00 - O*aO*Do10
# dDo01  = -KO*aO*Do01 + K_nt*a_nt*Do11 - NT*a_nt*Do01 + O*aO*Do00
# dDo11  = -KO*aO*Do11 - K_nt*a_nt*Do11 + NT*a_nt*Do01 + O*aO*Do10
# dO     = m3 + Do11*alphaO - O*delta + KO*aO*Dn1 + KO*aO*Do01 + KO*aO*Do11 + KO*aO*Dt01 + KO*aO*Dt11 - O*aO*Dn0 - O*aO*Do00 - O*aO*Do10 - O*aO*Dt00 - O*aO*Dt10
# dDt00  = KO*aO*Dt01 + K_nt*a_nt*Dt10 - NT*a_nt*Dt00 - O*aO*Dt00
# dDt10  = KO*aO*Dt11 - K_nt*a_nt*Dt10 + NT*a_nt*Dt00 - O*aO*Dt10
# dDt01  = -KO*aO*Dt01 + K_nt*a_nt*Dt11 - NT*a_nt*Dt01 + O*aO*Dt00
# dDt11  = -KO*aO*Dt11 - K_nt*a_nt*Dt11 + NT*a_nt*Dt01 + O*aO*Dt10
# dDn0   = KO*aO*Dn1 - O*aO*Dn0
# dDn1   = -KO*aO*Dn1 + O*aO*Dn0
# dDm    = Do00*gamma - NT*Dm*theta
