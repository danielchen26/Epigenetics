using DifferentialEquations, DiffEqBiological
using Plots;pyplot()

# Making Animation
"""
    anim_gen(u0, p, tspan, model, diff, max, idx, vars, ylim_max, t_max_plt, mi = 1)
...
# Parameters:
- `u0`                  : Inital condition.
- `p`                   : Model Parameter
- `model`               : Model selected
- `diff`                : The increase unit of initial variables Concentration
- `max`                 : The max values of initial variables Concentration
- `idx`                 : The index array indicates which variable will be initial turned on with non zero values
- `idx_o`               : The index indicates which variable should be over express.
- `OE`                  : The inital value of the over expressed variable.
- `vars`                : The final variables to be displayed
- `ylim_max`            : plot ylim max value
- `t_max_plt`           : t_max for plot
- `mi`                  : maxiters
...
"""
function anim_gen(u0,p, ts, model, diff, max, idx, idx_o, OE, vars,  mi = 1e6) # ylim_max, t_max_plt,

    anim = @animate for i=0:diff:max
        u0[idx] .= i # [3,4,10,13,14]
        u0[idx_o] .= OE
        oprob = ODEProblem(model, u0, ts, p) # O_N_auto
        osol  = solve(oprob, Tsit5(), maxiters = mi)
        plot(osol, vars = getindex(model.syms, vars), legend = true, lw =1)
        # plot(osol, vars = getindex(model.syms, vars), tspan = (0,t_max_plt), ylims = (0,ylim_max),legend = true, lw =1)

        title!("Model Evolution : Strength Response : DM_initial = $i ")
        xaxis!("Time")
        yaxis!("Concentration")
    end
    return anim
end
