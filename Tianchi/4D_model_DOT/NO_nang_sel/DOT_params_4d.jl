using DiffEqBiological, DifferentialEquations
using Plots;gr()
using LinearAlgebra, DataFrames, Queryverse
include(pwd()*"/functions.jl")

# ================ CRN of demethylation on top of TF binding =============
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
    # Nanog(no self regulation)
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
@add_constraints Demethy_TF_MC begin
  Do00 + Do01 + Do10 + Do11 + D5mc + D5hmc  = 1
  Dt00 + Dt01 + Dt10 + Dt11 = 1
  Dn0 + Dn1 = 1
end

p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 1, 1, 1,  0., 0.05, 0.]
ss = steady_states(Demethy_TF_MC,p)
sort!(ss, by = x -> x[1]) # sort by Nanog
sb = stability(ss,Demethy_TF_MC,p)
sb2 = stability_tianchi(ss,Demethy_TF_MC,p,3)




# show DataFrame
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
    #
    # eigs = [JE_stability(i,rn1d,p)[2] for i in ss]
    # println("The eigs are:")
    # # @show eigs
    # @show hcat(eigs...)
end




# ------- Directly distant 3d [full CRN projection] ------------
kh = 0:0.1:10.   # effective θ1
beta = 0:0.1:10. # effective θ2
a_dn = 0:0.1:10. # effective γ

range = 0:0.1:20.
function DOT_params( KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta; range = 0:0.1:10., param = "a_dn")
    C_3d =similar(range)
    for i = eachindex(range)

        if param == "a_dn"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, range[i], kh, beta, 0., 0.05, 0.]
        elseif param  == "kh"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, range[i], beta, 0., 0.05, 0.]
        elseif param == "beta"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, range[i], 0., 0.05, 0.]
        end

        @show range[i]
        ss = steady_states(Demethy_TF_MC,p)
        sort!(ss, by = x -> x[1])

        if length(ss) >2
            NTO_idx = [1,2,9]
            Low = ss[1][NTO_idx]; Mid = ss[2][NTO_idx]; High = ss[3][NTO_idx]
            @show Low Mid High

            # C_3d[i] = DOT = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3]))
            C_3d[i] = DOT = (norm(ss[1][NTO_idx]-ss[2][NTO_idx]))/(norm(ss[1][NTO_idx]-ss[3][NTO_idx]))
            @show DOT
        else
            C_3d[i] = NaN
        end
    end
    return C_3d
end


# function Plot_DOT(p, p_range; p_name = "a_dn")
#     C_3d = DOT_params(p..., range = p_range, param = p_name)
#     plot(p_range,C_3d, xlabel = p_name, ylabel = "Domain of Attraction %", legend = false)
# end
# p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1., 1., 1., 1.]
# Plot_DOT(p, range, p_name = "beta")


#  plot for parameter set
using ProgressMeter
pp = plot()
@showprogress for beta = 1:10
    p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1., 1., 1., beta]
    C_3d = DOT_params(p..., range = range, param = "a_dn")
    plot!(pp, range, C_3d, xlabel = "a_dn", ylabel = "Domain of Attraction %", legend = false)
end

pp








# Define the reduced ODE `Demethy_TF_ode_reduced_MC`  for the sampling below # Cons .= 1
using ParameterizedFunctions
Demethy_TF_ode_reduced_MC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O
    dN = m1 - N*delta + (O*alphaN)/(KO + O)
    dT = m2 - T*delta + (N^2*O*T^2*alphaT)/((2*K_nt*Kd^2 + N^2*T^2)*(KO + O))
    dO = (N^3*O*T^3*alphaO*beta*kh + KO*N^3*T^3*beta*kh*m3 + N^3*O*T^3*beta*kh*m3 - N^3*O^2*T^3*beta*delta*kh + 2*KO*K_nt*Kd^3*a_dn*beta*m3 - 2*KO*K_nt*Kd^3*O*a_dn*beta*delta - KO*N^3*O*T^3*beta*delta*kh + 2*KO*K_nt*Kd^2*N*T*a_dn*kh*m3 + 2*KO*K_nt*Kd^2*N*T*beta*kh*m3 + 2*K_nt*Kd^2*N*O*T*beta*kh*m3 - 2*K_nt*Kd^2*N*O^2*T*beta*delta*kh - 2*KO*K_nt*Kd^2*N*O*T*a_dn*delta*kh - 2*KO*K_nt*Kd^2*N*O*T*beta*delta*kh)/(2*KO*K_nt*Kd^3*a_dn*beta + KO*N^3*T^3*beta*kh + N^3*O*T^3*beta*kh + 2*KO*K_nt*Kd^2*N*T*a_dn*kh + 2*KO*K_nt*Kd^2*N*T*beta*kh + 2*K_nt*Kd^2*N*O*T*beta*kh)
end KO K_nt Kd a1 d a_nt aO alphaT alphaO alphaN delta a_dn beta kh m1 m2 m3

range = 0.:1.:10.
# ---------- Sample from 3d space (NTO) using Demethy_TF_ode_reduced_MC --------
function DOT_Volume_params(KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta; range = 0:0.1:10., param = "a_dn")
    C = similar(range)
    for i in eachindex(range)# add @suppress
        if param == "a_dn"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, range[i], kh, beta, 0., 0.05, 0.]
        elseif param  == "kh"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, range[i], beta, 0., 0.05, 0.]
        elseif param == "beta"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, range[i], 0., 0.05, 0.]
        end
        # Using the full CRN to check the # of SSS
        ss = steady_states(Demethy_TF_MC,p)
        sort!(ss, by = x -> x[1])

        if length(ss) >2
            NTO_idx = [1,2,9]
            Low = ss[1][NTO_idx]; Mid = ss[2][NTO_idx]; High = ss[3][NTO_idx]
            # @show Low Mid High

            smpl_max = extrema(vcat([each[NTO_idx] for each in ss]...))[2]*1.5
            cube_NTO = LinRange(0.,smpl_max,10)
            tspan = (0., 5e2)
            soma = []
            TV = 0
            for  N = cube_NTO, T = cube_NTO ,O = cube_NTO
                u0 = [N,T,O]
                # @show u0
                prob = ODEProblem(Demethy_TF_ode_reduced_MC,u0,tspan,p)
                sol = solve(prob,Rosenbrock23())
                f_ss = norm(sol[end] .- Low) < 0.1 ? 1 : 0
                # @show f_ss
                push!(soma, f_ss)
                TV += 1
                # @show sol[end]
            end
            # for i = 1:2e3
            #     u0 = init_rand(cube_O,con)
            #     prob = ODEProblem(rn1d_leak,u0,tspan,p)
            #     sol = solve(prob,Rosenbrock23())
            #     @show sol[end]
            #     f_ss = norm(sol[end] .- ss[1]) < 0.1 ? 1 : 0
            #     @show f_ss
            #     push!(soma, f_ss)
            # end
            DOT = sum(soma)/TV
            C[i] = DOT
        else
            C[i] = NaN
        end
    end
    return C
end


p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1., 1., 1., 1.]
p = [0.3, 0.2, 0.1, 1., 1., 10., 10., 1.0, 1.0, 1.0, 1., 1., .3, .3]
C = DOT_Volume_params(p..., range = [0:0.01:0.1; .1:.1:6.], param = "kh")
plot([0:0.01:0.1; .1:.1:6.], C)




# ========= DOT_Volume_params plot ===========
pv = plot()
range = 0:10.
@showprogress for KO = .3, K_nt = 0.2, Kd = 0.1, a1 = 1:1.:2.0, d = 1:1.:2.0, a_nt = 10, aO = 10, alphaT = 0:2., alphaO = 0:2., alphaN = 0:2., a_dn = 1., delta = 0:2., kh = 0:2., beta= 0:2.
    p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta]
    C = DOT_Volume_params(p..., range = range, param = "a_dn")
    @show C
    if sum(isnan.(C)) > 1
        continue
    else
        plot!(pv,range,C)
    end

    sleep(0.01)
end
plot(pv,xlabel="a_dn", ylabel = "Domain of Attraction %", legend = false)
savefig(pv,"/Users/chentianchi/Desktop/αβϵθδ.png")





pv = plot()
range = [0:0.01:0.1; .1:.1:10.]
@showprogress for KO = .3, K_nt = 0.2, Kd = 0.1, a1 = 1., d = 1., a_nt = 10., aO = 10., alphaT = 1., alphaO = 1., alphaN = 1.,  delta = 1., a_dn = 0:0.1:1., kh = 1., beta= 0:0.1:1.
    p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta]
    C = DOT_Volume_params(p..., range = range, param = "beta")
    @show C
    if sum(isnan.(C)) < length(C)
        plot!(pv,range,C, legend = false)
    else
        continue
    end

    sleep(0.01)
end

plot(pv,xlabel="beta", ylabel = "Domain of Attraction %", legend = false)
savefig(pv,"/Users/chentianchi/Desktop/full_a_dn_kh.png")







# ==== A 3d polygon separatrix search =======
function separatrix( KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, beta; range = 0:0.1:10., param = "a_dn")
    C = similar(range)
    for i in eachindex(range)# add @suppress
        if param == "a_dn"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, range[i], kh, beta, 0., 0.05, 0.]
        elseif param  == "kh"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, range[i], beta, 0., 0.05, 0.]
        elseif param == "beta"
            p = [KO, K_nt, Kd, a1, d, a_nt, aO, alphaT, alphaO, alphaN, delta, a_dn, kh, range[i], 0., 0.05, 0.]
        end
        # Using the full CRN to check the # of SSS
        ss = steady_states(Demethy_TF_MC,p)
        sort!(ss, by = x -> x[1])
        if length(ss) >2
            NTO_idx = [1,2,9]
            Low = ss[1][NTO_idx]; Mid = ss[2][NTO_idx]; High = ss[3][NTO_idx]
            # @show Low Mid High
            saddle = Mid
            smpl_max = extrema(vcat([each[NTO_idx] for each in ss]...))[2]*1.5 # model range
            cube = 0.1
            tspan = (0., 5e2)
            # Test the cube 6 points to
            for  N = cube_NTO, T = cube_NTO ,O = cube_NTO
                u0 = [N,T,O]
                # @show u0
                prob = ODEProblem(Demethy_TF_ode_reduced_MC,u0,tspan,p)
                sol = solve(prob,Rosenbrock23())
                f_ss = norm(sol[end] .- Low) < 0.1 ? 1 : 0
                # @show f_ss
                push!(soma, f_ss)
                TV += 1
                # @show sol[end]
            end
            # for i = 1:2e3
            #     u0 = init_rand(cube_O,con)
            #     prob = ODEProblem(rn1d_leak,u0,tspan,p)
            #     sol = solve(prob,Rosenbrock23())
            #     @show sol[end]
            #     f_ss = norm(sol[end] .- ss[1]) < 0.1 ? 1 : 0
            #     @show f_ss
            #     push!(soma, f_ss)
            # end
            DOT = sum(soma)/TV
            C[i] = DOT
        else
            C[i] = NaN
        end
    end
    return C
end











# ==== test time scale to equilibrate =======
# cube_NTO = LinRange(0.,15.,5)
# tspan = (0., 5e1)
# p = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1., 1., 1., 1.,   0., 0.05, 0.]
# anim = @animate for  N = cube_NTO, T = cube_NTO ,O = cube_NTO
#     u0 = [N,T,O]
#     @show u0
#     prob = ODEProblem(Demethy_TF_ode_reduced_MC,u0,tspan,p)
#     sol = solve(prob,Rosenbrock23())
#     plot(sol, vars = [:N, :T, :O], ylim = (0.,5.))
# end
# gif(anim, "/tmp/anim_fps15.gif", fps = 3)
