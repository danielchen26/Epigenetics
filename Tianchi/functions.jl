# ======= Collection of functions ========
using DifferentialEquations
# Random numbers with fixed sum function
function randfixsum(n, dim, tot)
    sol = []
    for i in 1: n
        v = rand(1,dim)
        nv = tot*v./sum(v)
        push!(sol, nv)
    end
    return sol
end


# function make_cb2(var1, var1_val, var2, var2_val,ts_in)
#     ts = ts_in
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#       if integrator.t == ts[1]
#           integrator.p[var1] = var1_val
#           integrator.p[var2] = var2_val
#       elseif integrator.t == ts[2]
#           integrator.p[var1] = 0.0
#           integrator.p[var2] = 0.0
#       end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     return ts, cb
# end


# function make_cb2(ts_in, vars...)
#     ts = ts_in
#     condition(u,t,integrator) = t in ts
#     function affect!(integrator)
#       if integrator.t == ts[1]
#           integrator.p[vars[1]] = vars[2]
#           integrator.p[vars[3]] = vars[4]
#           integrator.p[vars[5]] = vars[6]
#       elseif integrator.t == ts[2]
#           integrator.p[vars[1]] = 0.0
#           integrator.p[vars[3]] = 0.0
#           integrator.p[vars[5]] = 0.0
#       end
#     end
#     cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
#     @show vars[1], vars[2], vars[3], vars[4], vars[5], vars[6]
#     return ts, cb
# end






function make_cb(ts_in, vars...)
    ts = ts_in
    condition(u,t,integrator) = t in ts
    function affect!(integrator)
        for  i = 1:2: length(vars)
            if integrator.t == ts[1]
                integrator.p[vars[i]] = vars[i+1]
            elseif integrator.t == ts[2]
                integrator.p[vars[i]] = 0.0
            end
        end
    end
    cb = DiscreteCallback(condition, affect!, save_positions=(true,true));
    @show vars
    return ts, cb
end




function JE_stability(solution::Vector{Float64}, rn::DiffEqBase.AbstractReactionNetwork, p::Vector{Float64}, t=0.::Float64)
    jac = zeros(length(rn.syms),length(rn.syms))
    rn.jac(jac,solution,p,t)
    return (jac,eigen(jac).values)
end


using Reduce
function subs(old, new, Eq)
    Algebra.sub(old => new, Eq) |> mat
end
function solve(Eq, var)
    Algebra.solve(:(0  == $Eq),var)
end
