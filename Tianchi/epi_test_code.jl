using DifferentialEquations, ParameterizedFunctions
using Plots
using Latexify

module epi
    include("functions.jl")
    include("./Model_collection/ODE_model_set.jl")
    include("./Old_Unused/Init.jl")
end




# 1. ============ Single CRN test =============
# randomize Num initial condition for same conservation class
Num = 200
DO_tot = epi.randfixsum(Num,4)
DT_tot = epi.randfixsum(Num,4)
# ------ CRN initial condition ---------
@show u0, p = epi.crn_init(DO_tot,DT_tot,2)
tspan = (0.0,1.5e3)
# ------ Control(choose one) -------
# 1. Oct4 overexpress
ts, cb = epi.make_cb([500,7e2],13,0.06,15,0.36)
# 2. TET overexpress
ts, cb = epi.make_cb([500,7e2],13,0.06,14,0.36,15,0.06) # need modification
# 3. Nanog overexpress
ts, cb = epi.make_cb([500,7e2],13,0.36,15,0.06)
# ------ Solving
prob = ODEProblem(epi.paper_ode,u0,tspan,p)
sol = solve(prob, Rosenbrock23(),callback=cb, tstops=ts)
plot(sol, vars = [:N, :T, :O],  title = "CRN Model", ylabel= "Concentration", lw = 1.5, ylims = (0,1.5))







# 2. ============ Randomized initial condition CRN test ============
OTN_SS = []； OTN_SS_all = []
anim = @animate for i  = 1:Num
    # ------- initial condition -----
    u0, p = epi.crn_init(DO_tot,DT_tot,i)
    tspan = (0.0,1e4)

    # ------ Control(choose one) -------
    # 1. Oct4 overexpress
    # ts, cb = epi.make_cb([5e3,7e3],13,0.06,15,0.36)
    # 2. TET overexpress
    ts, cb = epi.make_cb([5e3,7e3],13,0.06,14,0.36,15,0.06)
    # 3. Nanog overexpress
    # ts, cb = epi.make_cb([5e3,7e3],13,0.36,15,0.06)

    # ------ Solving -------
    prob = ODEProblem(epi.paper_ode,u0,tspan,p)
    sol = solve(prob, DynamicSS(Rodas5()),callback=cb, tstops=ts)
    push!(OTN_SS, sol[end][[1,2,9]]); push!(OTN_SS_all, sol[end])
    # ------ Plotting ------
    plot(sol, vars = [:N, :T, :O],  w = 1.5, ylims = (0,1.5),
    title = "CRN with random initial $i", xlabel = "Time", ylabel= "Concentration")
end
gif(anim, "/Users/chentianchi/Desktop/Oct4_random_init.gif", fps = 5)







# ===================================   test reduced model
using DifferentialEquations, ParameterizedFunctions, LinearAlgebra, Plots, Measures
module epi
    include("functions.jl")
    include("ODE_model_set.jl")
    include("Init.jl")
end
gr(dpi =200)

paper_reduce_ode = @ode_def begin
    dN = -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1
    dT = (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    dO = (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3


p = [0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6]
tspan = (0.0,1e4)
ts, cb = epi.make_cb([5e3,7e3],13,0.06,15,0.36) #oct4_over
prob_rd = ODEProblem(paper_reduce_ode,[0.,0.1,0.],tspan,p)
sol_rd = solve(prob_rd, DynamicSS(Rodas5()),callback=cb, tstops=ts)
sol_rd = solve(prob_rd, Rosenbrock23(),callback=cb, tstops=ts, abstol=1e-8)
plot(sol_rd, w = 1.5, ylims = (0,1.5), title = "Reduced Model N T O", xlabel = "Time", ylabel= "Concentration")

# anim
anim = @animate for i = 1:100
    prob_rd = ODEProblem(paper_reduce_ode,rand(3),tspan,p)
    sol_rd = solve(prob_rd, DynamicSS(Rodas5()),callback=cb, tstops=ts, abstol=1e-8)
    plot(sol_rd, w = 1.5, ylims = (0,1.5), title = "Reduced Model N T O", xlabel = "Time", ylabel= "Concentration")
end
gif(anim, "/tmp/rd_NTO.gif", fps = 5)




#  -------   Parameter explortion
# Preset values Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3
p = [0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6]

tspan = (0.0,5e2)
ts, cb = epi.make_cb([3e2,4e2],13,0.06,15,0.36) #oct4_over

# -------parameters with random initial condition
anim = @animate for p[1] = 0.3:0.05:0.4, p[2] = 0.2:0.05:0.3, p[3] = 0.1:0.05:0.2
    @show p
    # create a plot with subplots and a custom layout
    plot_array = Any[]
    for i  = 1:6
        u0 = rand(3)
        prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
        sol_rd = solve(prob_rd, Rosenbrock23(),callback=cb, tstops=ts, abstol=1e-8)

        t = 0.0
        out = rand(3,3)
        paper_reduce_ode.jac(out,sol_rd[end],p,t)
        @show out
        eig = eigen(out).values
        push!(plot_array,bar(eig, legend = false))

        push!(plot_array,plot(sol_rd, w = 1.5, xlabel= " ", ylims = (0,1.5),  title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 7, tickfontsize = 4, legend = false))
    end
    plot(plot_array...,layout = (6,2), margin =0mm)
end
gif(anim, "/tmp/K_dissociation_V.gif", fps = 1)




#  -------- save for each parameter set
anim = @animate for p[1] = 0.3:0.05:0.4, p[2] = 0.2:0.05:0.3, p[3] = 0.1:0.05:0.2
    # @show p
    # create a plot with subplots and a custom layout
    plot_array = Any[]
    for i  = 1:35
        @show u0 = rand(3)
        prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
        sol_rd = solve(prob_rd, Rosenbrock23(),callback=cb, tstops=ts, abstol=1e-8)

        t = 0.0
        out = rand(3,3)
        paper_reduce_ode.jac(out,sol_rd[end],p,t)
        @show out
        eig = eigen(out).values
        @show eig
        push!(plot_array,bar(eig, title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4, tickfontsize = 4,legend = false))

        # push!(plot_array,plot(sol_rd, w = 1.5, xlabel= " ", ylims = (0,1.5),  title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 7, tickfontsize = 4, legend = false))
    end
    plot(plot_array...,layout = (7,5), margin =0mm)
end
gif(anim, "/tmp/K_dissociation_V.gif", fps = 1)





# same random initial condition with different set of parameters
# animation
anim = @animate for i = 1:50
    u0 = rand(3)
    plot_array_ds = Any[]
    # plot_array_eig = Any[]
    out = rand(3,3)
    for p[1] = 0.3:0.05:0.4, p[2] = 0.2:0.05:0.3, p[3] = 0.1:0.05:0.2
        # @show p

        prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
        sol_rd = solve(prob_rd, Rosenbrock23(),callback=cb, tstops=ts, abstol=1e-8)

        t = 0.0
        paper_reduce_ode.jac(out,sol_rd[end],p,t)
        # @show out
        eig = eigen(out).values
        @show eig
        # push!(plot_array_eig,bar(eig, title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4, tickfontsize = 4,legend = false, c = :viridis, clims = [-1,0]))

        push!(plot_array_ds,plot(sol_rd, w = 1.5,  xlabel= " ", ylims = (0,1.5),  title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4,  tickfontsize = 4, legend = false))
    end
    # plot(plot_array_eig ...,  margin = 0mm)
    plot(plot_array_ds...,  margin = 0mm)
end
gif(anim, "/Users/chentianchi/Desktop/param_C.gif", fps = 2)



# save each DS plots with same_init
#  -----  DS plot
tspan = (0.0,5e2)
ts, cb = epi.make_cb([3e2,4e2],13,0.06,15,0.36) #oct4_over
for i = 1:20
    u0 = rand(3)
    plot_array_ds = Any[]
    # plot_array_eig = Any[]
    for p[1] = 0.3:0.05:0.4, p[2] = 0.2:0.05:0.3, p[3] = 0.1:0.05:0.2
        @show p
        prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
        sol_rd = solve(prob_rd, Rosenbrock23(),callback=cb, tstops=ts, abstol=1e-8)

        t = 0.0
        out = rand(3,3)
        paper_reduce_ode.jac(out,sol_rd[end],p,t)
        @show out
        eig = eigen(out).values
        # push!(plot_array_eig,bar(eig, title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4, tickfontsize = 4,legend = false, c = :viridis, clims = [-1,0]))

        push!(plot_array_ds,plot(sol_rd, w = 1.5, xlabel= " ", ylims = (0,1.5),  title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4, tickfontsize = 4, legend= false))
    end
    plt_ds = plot(plot_array_ds..., margin =0mm)
    # plt_eig = plot(plot_array_eig..., margin =0mm)

    savefig(plt_ds, "/Users/chentianchi/Desktop/epi_plots/params_explore/julia/same_init/DS_$i.png")
    # savefig(plt_eig, "/Users/chentianchi/Desktop/epi_plots/params_explore/julia/same_init/Eig_$i.png")
end





# -------- save each Eig plot for same init
tspan = (0.0,1e3)
ts, cb = epi.make_cb([3e2,4e2],13,0.06,15,0.36) #oct4_over
for i = 1:20
    u0 = rand(3)
    # plot_array_ds = Any[]
    plot_array_eig = Any[]
    for p[1] = 0.3:0.05:0.4, p[2] = 0.2:0.05:0.3, p[3] = 0.1:0.05:0.2
        @show p
        prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
        sol_rd = solve(prob_rd, Rosenbrock23(),callback=cb, tstops=ts, abstol=1e-8)

        t = 0.0
        out = rand(3,3)
        paper_reduce_ode.jac(out,sol_rd[end],p,t)
        @show out
        eig = eigen(out).values
        push!(plot_array_eig,bar(eig, title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4, tickfontsize = 4,legend = false, c = :viridis, clims = [-1,0]))

        # push!(plot_array_ds,plot(sol_rd, w = 1.5, xlabel= " ",  ylims = (0,1.5),  title = "K_o : $(p[1])  K_nt : $(p[2])   Kd : $(p[3])", titlefontsize = 4, tickfontsize = 4, legend= false))
    end
    # plt_ds = plot(plot_array_ds..., margin =0mm)
    plt_eig = plot(plot_array_eig..., margin =0mm)

    # savefig(plt_ds, "/Users/chentianchi/Desktop/epi_plots/params_explore/julia/same_init/DS_$i.png")
    savefig(plt_eig, "/Users/chentianchi/Desktop/epi_plots/params_explore/julia/same_init/Eig_$i.png")
end




# ================ eig_spectrum ==============
p = [0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, .05, 1e-6]
tspan = (0.0,1.5e2)
ts, cb = epi.make_cb([5e1,7e1],13,0.06,15,0.36) #oct4_over

eig_spectrum = []
for i  = 1:50
    u0 = rand(3)
    prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(), tstops=ts, abstol=1e-8)

    t = 0.0
    out = rand(3,3)
    paper_reduce_ode.jac(out,sol_rd[end],p,t)
    display(out)
    eig = eigen(out).values
    display(eig)
    b_eig = bar(eig,legend = false, tickfontsize = 2,c = :viridis, clims = [-1,0])
    push!(eig_spectrum, b_eig)
end

plot(eig_spectrum...)























# --- Sampling the space ----- trying a volume plot in 3D space.
using DifferentialEquations, ParameterizedFunctions, LinearAlgebra, Measures, Queryverse
using DataFrames, Makie, Images #Plots, VegaLite
using Distributions, ProgressMeter
module epi
    include("functions.jl")
    include("./Model_collection/ODE_model_set.jl")
end


paper_reduce_ode = @ode_def begin
    dN = -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1
    dT = (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    dO = (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
end Kₒ Kₙₜ Kd a₁ d aₙₜ aₒ αₜ αₒ αₙ δ β m1 m2 m3


p = [0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, .05, 1e-6]
u0 = [0.304046, 0.181063, 0.131063]; tspan = (0,1500.)
prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
@show sol_rd[end]
sol_rd(12)

# ----- for a 3D unit cubit, show SSS --------
# SSS_1 :=> [0.6953744055154728, 0.734880954107109, 0.6849131746074149]
# SSS_2 :=> [4.333320581187469e-6, 0.049999999998359725, 1.0000000000381673e-6]
N = 1; sp = 0.1; r = 0:sp:N
for x = r ,y =r, z=r
    u0 = [x,y,z]
    prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
    # @show sol_rd[end]
    @show sol_rd[end]
end

# Point_T_scaler calculate time it takes for each point in 3D fall into attractor
function Point_T_scaler(x,y,z; tspan = (0.0,1.5e2), p = p)
    u0 = [x,y,z]
    prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
    # @show sol_rd[end]
    t_lmin = sol_rd.t[findlocalminima(sol_rd[2,:])[end][1]]
    t_lmax = sol_rd.t[findlocalmaxima(sol_rd[2,:])[end][1]]
    t_ss = min(t_lmin, t_lmax)
    # @show t_ss
    sss = sol_rd[end]
    # @show sss
    attractor_high = abs(sum(sss .- 1)) < abs(sum(sss .- 0)) ? t_ss : 0
    attractor_low = abs(sum(sss .- 1)) > abs(sum(sss .- 0)) ? t_ss : 0
    return attractor_high, attractor_low
end



#  ------ show which attractor each point initial end ups to ---------
r = 0:0.5:1
@showprogress "Computing..." for x = r ,y =r, z=r
    attractor_high, attractor_low = Point_T_scaler(x,y,z)
    @show attractor_high
    @show attractor_low
    sleep(0.01)
end



r = 0:0.05:1
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])
high_select = (Tc_high .!= 0); low_select = (Tc_low .!= 0)
attract_high = @. high_select * (1/(Tc_high )); attract_low = @. low_select * (1/(Tc_low ))
hi_min, hi_max = sort(unique(attract_high))[[2,end-1]]
lo_min, lo_max = sort(unique(attract_low))[[2,end]]
p_high = volume(r,r,r, attract_high, algorithm = :mip, colorrange = (hi_min, hi_max))
p_low = volume(r,r,r, attract_low, algorithm = :mip, colorrange = (lo_min, lo_max))
cm_high = colorlegend(
           p_high[end],             # access the plot of Scene p1
           raw = true,          # without axes or grid
           camera = campixel!,  # gives a concrete bounding box in pixels
                                # so that the `vbox` gives you the right size
           width = (            # make the colorlegend longer so it looks nicer
               30,              # the width
               540              # the height
           )
           )
cm_low = colorlegend(
           p_low[end],             # access the plot of Scene p1
           raw = true,          # without axes or grid
           camera = campixel!,  # gives a concrete bounding box in pixels
                                # so that the `vbox` gives you the right size
           width = (            # make the colorlegend longer so it looks nicer
               30,              # the width
               540              # the height
           )
           )
scene_final = vbox(p_high,cm_high,p_low,cm_low)




# contour(r,r,r, Tc, alpha = 0.1, transparency = true, algorithm = :mip)








# heatmap slice control
s, value = textslider(1:size(Tc, 3), "slice"
hmap = heatmap(lift(idx-> Tc[:, :, idx], value), padding = (0.0, 0.0))
hbox(s, hmap)

# volume control
s, value = textslider(LinRange(extrema(Tc)..., 100), "isovalue")
density = Makie.volume(Tc,isovalue = value,  algorithm = :iso)
hbox(s, density)







# ---- ugly code
# function f_mat(tspan,p,N,sp)
#     mat = []
#     for i= 0:sp:N, j= 0:sp:N, k= 0:sp:N
#         # @show i j k
#         u0 = [i,j,k]
#         prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
#         sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
#         # @show sol_rd[end]
#         t_lmin = sol_rd.t[findlocalminima(sol_rd[2,:])[end][1]]
#         t_lmax = sol_rd.t[findlocalmaxima(sol_rd[2,:])[end][1]]
#         t_ss = min(t_lmin, t_lmax)
#         @show t_ss
#         sss = sol_rd[end]
#         id = abs(sum(sss .- 1)) < abs(sum(sss .- 0)) ? t_ss : 0
#         push!(mat, i,j,k,id)
#     end
#     return data = DataFrame(reshape(mat, 4, Int((N/sp + 1)^3))')
# end
# # make above a interpolatable function??-------
# data =f_mat(tspan,p,N,sp)
#
# function f_match(data, x, y, z)
#     scaler_t = @from i in data begin
#         @where  i.x1 == x && i.x2 ==y && i.x3 ==z
#         @select i.x4
#         @collect DataFrame
#     end
#     return Matrix(scaler_t)[2]
# end
