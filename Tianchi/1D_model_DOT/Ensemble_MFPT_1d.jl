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
    Di + Da + Dm = 50.
end


# using Parameters
mutable struct rn{T}
    α::T; ϵ::T
    β::T; η::T
    γ::T; θ::T
    δ::T
end

gr()
plotly()
using Blink
Plots.scalefontsizes(2)
w = Window()


function MFPT( model, p::rn; tspan = (0.,1e3))
    @show params = [p.α, p.ϵ, p.β, p.η, p.γ, p.θ, p.δ]

    ss = steady_states(model,params)
    sort!(ss, by = x -> x[1])

    if length(ss) >2
        @show ss[1];@show ss[2];@show ss[3]
        # start from metastable state
        u0_middle = floor.(ss[2])
        @show u0_middle
        @show conservation = sum(u0_middle[2:end])

        dprob_middle = DiscreteProblem(model, u0_middle, tspan, params)
        jprob_middle = JumpProblem(dprob_middle, Direct(), model)
        jsol_middle = solve(jprob_middle, SSAStepper())

        p1_middle = plot(jsol_middle, lw = 0.5, title = "SDE Start from Metastable");
        p2_middle = histogram(jsol_middle[1,:], bin = 50, legend = false, title = "Oct4 histogram");

        Cs_middle = map(x -> x < ss[2][1],jsol_middle[1,:] )
        p3_middle = histogram(jsol_middle[1,:][Cs_middle], legend = false, title = "Somatic Domain")
        Cl_middle = map(x -> x > ss[2][1],jsol_middle[1,:] )
        p4_middle = histogram(jsol_middle[1,:][Cl_middle], legend = false, title = "iPS Domain")

        l = @layout [a ; b  [c; d]]
        ui_middle = plot(p1_middle, p2_middle, p3_middle, p4_middle, layout = l)

        # start from iPS state
        u0_high = floor.(ss[3])
        @show u0_high
        @show conservation = sum(u0_high[2:end])

        dprob_high = DiscreteProblem(model, u0_high, tspan, params)
        jprob_high = JumpProblem(dprob_high, Direct(), model)
        jsol_high = solve(jprob_high, SSAStepper())

        p1_high = plot(jsol_high, lw = 0.5, title = "SDE Start from iPS");
        p2_high = histogram(jsol_high[1,:], bin = 50, legend = false, title = "Oct4 histogram");

        Cs_high = map(x -> x < ss[2][1],jsol_high[1,:] )
        p3_high = histogram(jsol_high[1,:][Cs_high], legend = false, title = "Somatic Domain")
        Cl_high = map(x -> x > ss[2][1],jsol_high[1,:] )
        p4_high = histogram(jsol_high[1,:][Cl_high], legend = false, title = "iPS Domain")


        l = @layout [a ; b  [c; d]]
        ui_high = plot(p1_high, p2_high, p3_high, p4_high, layout = l)

        l2 = @layout [m n]
        ui = plot(ui_middle, ui_high, layout = l2)


        body!(w, ui)
    else
        println("Single State")
    end
end


p =rn(0.5, 5., 5.,0.1,5.,0.2,5.)
MFPT(rn1d_leak, p::rn)

anim = @animate for γ = 4:0.1:5
    p =rn(0.5, 5., 5.,0.1,γ,0.2,5.)
    MFPT(rn1d_leak, p::rn)
end
gif(anim, "/tmp/test.gif", fps = 1)




u0 = [1., 2., 1., 10.]
tspan = (0.,1e3)
params = [5., 5., 5.,0.1,5.,5.,5.]
dprob = DiscreteProblem(rn1d_leak, u0, tspan, params)
jprob = JumpProblem(dprob, Direct(), rn1d_leak)
jsol = solve(jprob, SSAStepper())
ui= plot(jsol)




l = @layout [  a{0.3w} [grid(3,3)
                         b{0.5h} ]]
plot(
    rand(10,11),
    layout = l, legend = false, seriestype = [:bar :scatter ],
    title = ["($i)" for j = 1:1, i=1:11], titleloc = :right, titlefont = font(8)
)












# ====== just SSS display ====
@manipulate for α = 0:0.1:10.0 , ϵ = 0:0.1:10., β = 0:0.1:10.0, η  = 0:0.1:3.0, γ = 0:0.1:10.,  θ= 0:0.1:10., δ = 0:0.1:10.0
    p = [α, ϵ, β, η, γ, θ, δ]
    ss = steady_states(rn1d_leak,p)
    sort!(ss, by = x -> x[1])
    @show ss
end












using Interact;plotly()
@manipulate for α = 0:0.1:10.0 , ϵ = 0:0.1:10., β = 0:0.1:10.0, η  = 0:0.1:3.0, γ = 0:0.1:10.,  θ= 0:0.1:10., δ = 0:0.1:10.0
    p = [α, ϵ, β, η, γ, θ, δ]
    ss = steady_states(rn1d_leak,p)
    sort!(ss, by = x -> x[1])

    if length(ss) >2

        @show ss

        u0 = [1.,10.,0.,40.]
        tspan = (0., 1e3)
        dprob = DiscreteProblem(rn1d_leak, u0, tspan, p)
        jprob = JumpProblem(dprob, Direct(), rn1d_leak)
        jsol = solve(jprob, SSAStepper())
        plot(jsol) # vars=[:O, :Da]
    else
        @show "single state"
    end
end



α, ϵ, β, η, γ, θ, δ = 5.2 ,9.7, 9.9,0.1,5.8,0.3,9.4
p = [α, ϵ, β, η, γ, θ, δ]
u0 = [1.,10.,10.,30.]
tspan = (0., 1e3)
dprob = DiscreteProblem(rn1d_leak, u0, tspan, p)
jprob = JumpProblem(dprob, Direct(), rn1d_leak)
jsol = solve(jprob, SSAStepper())
plot(jsol)
df = hcat(jsol.u...)'
df[:,1]

scattermat(df[:,[1,3]], dims = 1)

cor2(df[:,1], df[:,3])

using Statistics
cor([rand(1000000,2);rand(1000000,2) .+ 1])
