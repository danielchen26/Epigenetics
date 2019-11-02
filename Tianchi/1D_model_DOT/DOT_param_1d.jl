using DiffEqBiological, DifferentialEquations
using Plots;gr()
using LinearAlgebra
module fun
    include("../functions.jl")
end


# 1. ======== CRN ========

rn1d_leak = @reaction_network begin
 (α, ϵ),            2O + Di ↔ Da
  β,                Da → Da + O
  η,                Di → Di + O
 (γ,θ),             Di ↔ Dm
  δ,                O → ∅
end α ϵ β η γ θ δ

@add_constraints rn1d_leak begin
    Di + Da + Dm = 5.
end
p = [1., 0.2, 1., 0.03, 1., 0.25, 0.22]
p = [5., 5., 5., 0.1, 10.5, 5., 5.]
ss = steady_states(rn1d_leak,p)
sort!(ss, by = x -> x[1])




# First Visulization =====================
# ===== 1d model DOT defined by Oct4 -----
using Interact
@manipulate for α = 0:0.1:10.0 , ϵ = 0:0.1:10., β = 0:0.1:10.0, η  = 0:0.1:3.0, γ = 0:0.1:10.,  θ= 0:0.1:10., δ = 0:0.1:10.0
    p = [α, ϵ, β, η, γ, θ, δ]
    ss = steady_states(rn1d_leak,p)
    sort!(ss, by = x -> x[1])
    if length(ss) >2
        DOT_proj = (norm(ss[1][1]-ss[2][1]))/(norm(ss[1][1]-ss[3][1]))
        DOT      = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3]))
        @show DOT_proj, DOT
    else
        @show "single state"
    end
    plot(sort([i[1] for i in ss]), marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
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





# ------- Directly distant 1d 4d ------------
using ProgressMeter,Suppressor
# ===== Define DOT -- γ --------
γ_rg = 0:0.1:10.
function DOT_γ( α, ϵ, β, η, θ, δ)
    C_1d =similar(γ_rg)
    C_3d =similar(γ_rg)
    for i = eachindex(γ_rg)
        p = [α, ϵ, β, η, γ_rg[i], θ, δ]
        ss = steady_states(rn1d_leak,p)
        sort!(ss, by = x -> x[1])
        if length(ss) >2
            C_1d[i] = DOT_proj = (norm(ss[1][1]-ss[2][1]))/(norm(ss[1][1]-ss[3][1]))
            C_3d[i] = DOT = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3]))
            # @show DOT_proj, DOT
        else
            C_1d[i] = NaN
            C_3d[i] = NaN
        end
    end
    return C_1d, C_3d
end




# ------- MC Sampling and DOT calculation --------- in 4d with Con = 5 --------
# # random initial functions
# function init_rand(cube_O,con)
#     u0 = [rand(cube_O)]
#     for i in fun.randfixsum(1,3,con)[1] push!(u0,i) end
#     return u0
# end
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


# ===test above two functions ====
C_1d, C_3d = DOT_γ(5., 5., 5., 0.1, 5., 5.)
C = DOT_Volume_params(5., 5., 5., 0.1, 5., 0., 5., range = [0:0.1:2.;2.:1.:10.], param = "θ")
plot([0:0.1:2.;2.:1.:10.],C)



#  ======= plot DOT_Volume_params sample ======
pp_sample = plot()
range = 0:20.
@showprogress for α = 0:0.5:5., ϵ = 0:0.5:5., β = 0:0.5:5., γ = 4.:0.5:5., δ = 4.:0.5:5.
    C = DOT_Volume_params(α, ϵ, β, 0.1, γ, 5., δ, range = range, param = "θ")
    if sum(isnan.(C)) < length(C)
        @show C
        plot!(pp_sample,range,C)
    else
        continue
    end
    sleep(0.01)
end

plot(pp_sample,xlabel="Theta", ylabel = "Domain of Attraction %", legend = false)
title!("DOT Sampling")
savefig(pp_sample,"/Users/chentianchi/Desktop/DOT_volume_θ_all.png")




# ========= DOT_γ plot ===========
pp = plot()
range = 0.:20.
@showprogress for α = 1.: 0.5: 5., ϵ = 1.: 0.5: 5., β = 1.: 0.5: 5.,θ = 4.:0.1:5.2, δ = 4.:0.1:5.2
    Ci_1d, Ci_3d = DOT_γ(α, ϵ, β, 0.1, θ, δ)
    # @show Ci_1d
    ind = Ci_1d[Ci_1d.>=0]
    plot!(pp,γ_rg[1:length(ind)],Ci_3d[1:length(ind)])
    sleep(0.01)
end
plot(pp,xlabel="Gama", ylabel = "Domain of Attraction %", legend = false)
savefig(pp,"/Users/chentianchi/Desktop/αβϵθδ.png")
# # ----- try to test some Parameters range
# using ProgressMeter, Suppressor
# P_set = []
# @showprogress for α = 0.:5.0 , ϵ = 0.:5.0, β = 0.:5.0, η  = 0.:3.0, γ = 0.:3.,  θ= 0.:2.0, δ = 0.:2.0
#     p = [α, ϵ, β, η, γ, θ, δ]
#     ss = steady_states(rn1d_leak,p)
#     if length(ss) >2
#         push!(P_set,p)
#     else
#         continue
#     end
#     sleep(0.01)
# end
# hcat(P_set...)












































# # ode ensemble_prob plots
# using Distributed
# using DifferentialEquations
# using Plots
#
# addprocs()
# @everywhere using DifferentialEquations
#
# # Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
# prob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))
#
# @everywhere function prob_func(prob,i,repeat)
#   ODEProblem(prob.f,rand()*prob.u0,prob.tspan)
# end
#
# ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
# sim = solve(ensemble_prob,Tsit5(),trajectories=100)
#
# plot(sim,la=0.9)


# # Solving an SDE with Different Parameters
# function f(du,u,p,t)
#   du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
#   du[2] = -3 * u[2] + u[1]*u[2]
# end
# function g(du,u,p,t)
#   du[1] = p[3]*u[1]
#   du[2] = p[4]*u[2]
# end
# p = [1.5,1.0,0.1,0.1]
# prob = SDEProblem(f,g,[1.0,1.0],(0.0,10.0),p)
#
# function prob_func(prob,i,repeat)
#   prob.p[3:4] = 0.3rand(2)
#   prob
# end
#
# ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
# sim = solve(ensemble_prob,SRIW1(),trajectories=10)
# plot(sim,linealpha=0.6,color=:blue,vars=(0,1),title="Phase Space Plot")
# plot!(sim,linealpha=0.6,color=:red,vars=(0,2),title="Phase Space Plot")
