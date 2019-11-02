# ================ CRN ======================
using DiffEqBiological, LinearAlgebra
using Plots;gr()
using DataFrames, Queryverse, Latexify

#
# # ----test 3 SSS -------
# simple_crn = @reaction_network begin
#     (ra,    ka*ra),   B2 + a0 ↔ a1
#     (ra,     ka*ra),  A2 + a0 ↔ a1
#     (rap, kap*rap),   a1 ↔ a1 + A
#     (rb,    kb*rb),   A2 + b0 ↔ b1
#     (rb,    kb*rb),   B2 + b0 ↔ b1
#     (rbp, kbp*rbp),   b1 ↔ b1 + B
#     (d,d),            B + B ↔ B2
#     (d,d),            A + A ↔ A2
#     (δ,δ),            (A, B) → ∅
# end ka ra kap rap kb rb kbp rbp d δ
#
# p = [1., 1., 10., 1.,1., 1., 10., 1., 10., .1]
# p = [5.,5.,5.,1.08,0.28,9.18,1.08,8.92,0.26,0.13]
# @add_constraints simple_crn begin
#   a0 + a1 = 1
#   b0 + b1 = 1
# end
#
# ss = steady_states(simple_crn,p)
# stability(ss,simple_crn,p)
#
# JE_stability(ss[2],simple_crn,p)[2]
#
#
#
#
#
#
#
#
#
#
# bif = bifurcation_grid_diagram(simple_crn, p,:kap,0.:10.,:kbp,(.1,5.))
# plot(bif)
# plot!([[],[]],color=[:blue :red],label = ["Stable" "Unstable"])
# # =============================
#
#
# sort!(ss, by = x -> x[4]) # sort by Nanog
# # show DataFrame
# dfc = DataFrame(vcat(ss))
# dfc.name = simple_crn.syms
# var = [:A, :B]
# @show dfc1 = dfc |> @filter(_.name in var) |> DataFrame




# First Visulization =====================
# ===== 2d model BOA defined by A-B-----

# -------- Jacobian -----------
function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,p,t)
    return (jac,eigen(jac).values)
end




# using Interact, Suppressor, ProgressMeter
# @manipulate for ka = 0:0.01:10.,ra = 0:0.01:10.,kap = 0:0.01:10.,rap = 0:0.01:10.,kb = 0:0.01:10.,rb = 0:0.01:10.,kbp = 0:0.01:10.,rbp = 0:0.01:10.,d = 0:0.01:10., δ =  0:0.01:10.
#     p = [ka, ra, kap, rap, kb, rb, kbp, rbp, d, δ]
#     ss = steady_states(simple_crn,p)
#     sort!(ss, by = x -> x[4])
#
#     dfc = DataFrame(vcat(ss))
#     dfc.name = simple_crn.syms
#     var = [:A, :B]
#     df_iPS = dfc |> @filter(_.name in var) |> DataFrame
#     @show Matrix(df_iPS)
#     # if length(ss) >2
#     #     DOT   = (norm(ss[1]-ss[2]))/(norm(ss[1]-ss[3]))
#     #     @show DOT
#     # else
#     #      @show "single state"
#     # end
#     eigs = [JE_stability(i,simple_crn,p)[2] for i in ss]
#     println("The eigs are:")
#     # @show eigs
#     @show hcat(eigs...)
#
#     plot(sort([i[[4,8]] for i in ss]), xlabel = "A-B",marker = (:hexagon, 10, 0.7, :green, stroke(1, 0.1, :black, :dot)))
#
# end


# adair(x1, x2,P::pp) = (P.a*x1^2 + P.b*x2^2 + P.c*x1^2*x2^2)/(1+ x1^2 + x2^2 + P.d*x1^2*x2^2)
# P1 = pp(0.276,  1.38, 0.897, 1.)
# P2 = pp(0.00828,0.0828, 0.092,1.)
# adair1(x1,x2) = adair(x1, x2,P1)
# adair2(x1,x2) = adair(x1, x2,P2)



# -------------- Below is collins paper model ---------------
# =======================# =======================# =======================
rn1 = @reaction_network begin
  (0.276*x1^2 + 1.38*x2^2 + 0.897*x1^2*x2^2)/(1+ x1^2 + x2^2 + x1^2*x2^2), ∅ → x1
  (0.00828*x1^2 + 0.0828*x2^2 + 0.092*x1^2*x2^2)/(1+ x1^2 + x2^2 + x1^2*x2^2), ∅ → x2
  # adair1(x1, x2), ∅ → x1
  # adair2(x1, x2), ∅ → x2
  d1, x1 → ∅
  d2, x2 → ∅
end d1 d2
p = [0.138, 0.046]
ss = steady_states(rn1,p)
stability(ss,rn1,p)


# --- if delete x1 self regulation ----
rn2 = @reaction_network begin
  (1.38*x2^2)/(1+  x2^2 ), ∅ → x1
  (0.00828*x1^2 + 0.0828*x2^2 + 0.092*x1^2*x2^2)/(1+ x1^2 + x2^2 + x1^2*x2^2), ∅ → x2
  # adair1(x1, x2), ∅ → x1
  # adair2(x1, x2), ∅ → x2
  d1, x1 → ∅
  d2, x2 → ∅
end d1 d2
p = [0.138, 0.046]
ss = steady_states(rn2,p)
stability(ss,rn2,p)
# =======================# =======================# =======================










# using ParameterizedFunctions
# struct pp{T}
#     a::T
#     b::T
#     c::T
#     d::T
# end
# Hill(x1, x2,P::pp) =  (P.a*x1^2 + P.b*x2^2 + P.c*x1^2*x2^2)/(1+ x1^2 + x2^2 + P.d*x1^2*x2^2)
#
# g1 = 0.138; g2 = 0.046
# P1 = pp(0.276,  1.38, 0.897, 1.)
# P2 = pp(0.00828,0.0828, 0.092,1.)
# DC_ode = @ode_def DC_ODE begin
#     dx1 = Hill(x1,x2,P1) - g1*x1
#     dx2 = Hill(x1,x2,P2) - g2*x2
# end P1 P2
# p = [P1,P2]
# u0 = [0.2,0.2]
# tspan = [0.,1e3]
# prob = ODEFunction(DC_ODE,u0,tspan,p)
#
# sol = solve(prob,SSRootfind())
# sol = solve(prob,DynamicSS(Tsit5()))
# using Sundials
# sol = solve(prob,DynamicSS(CVODE_BDF()),dt=1.0)
