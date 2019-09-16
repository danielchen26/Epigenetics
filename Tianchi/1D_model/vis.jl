
# ========================= Makie Vis ======================================
# ========================= Makie Vis ======================================
# ========================= Makie Vis ======================================

# ------------------------------------------------------------------------------------------------------------
# Eliminate Da ------------------

#  ----------- 3D vectorfield ----------
using Parameters, Makie, GeometryTypes
using LinearAlgebra
@with_kw struct rn1d_p{T}
    a::T
    ep::T
    b::T
    gama::T
    th::T
    d::T
end
P = rn1d_p(1., 0.2, 1., 0.25, 0.25, 0.22)
@unpack a, ep, b, gama, th, d = P



f(O,Di,Dm,P::rn1d_p) = Point3f0(
    -O*P.d + P.b*(1 -Di -Dm) + 2*P.ep*(1 -Di -Dm) - O^2*P.a*Di,
    -Di*P.gama + P.ep*(1 -Di -Dm) + P.th*Dm + (-1/2)*O^2*P.a*Di,
    Di*P.gama - P.th*Dm)

f(x) = f(x, P)


N = 20
r = LinRange(0., 5.5, N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))

vectorfield_dm = f.(x, y, z, Ref(P))
velocity_dm = norm.(vectorfield_dm)
normies_dm = sqrt.(velocity_dm) ./ 10.0
s1, a = textslider(0.07:0.1:0.27, "Green")
s2, b = textslider(0.07:0.1:0.27, "Yellow")
s3, c = textslider(0.07:0.1:0.27, "Red")
vfield_dm = meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield_dm),
    marker = Pyramid(Point3f0(0), 1f0, 1f0),
    markersize = Vec3f0.(0.2, 0.2, vec(normies_dm)),
    color = vec(normies_dm),
    # colormap = lift((x,y,z)-> [(:green, x), (:yellow, y), (:red, z)], a,b,c),
    colormap = [(:green, extrema(normies_dm)[1]), (:yellow, extrema(normies_dm)[2]), (:red, extrema(normies_dm)[2]+1)],
    colorrange = (extrema(normies_dm)[1], extrema(normies_dm)[2]),
    shading = false,
)
vbox(vfield_dm, colorlegend(vfield_dm[end]))
vbox(s1,s2,s3,vfield_dm, colorlegend(vfield_dm[end]))











# ------------- time to Equilibrate ------------
using DifferentialEquations, ParameterizedFunctions

rn1d_ode = @ode_def begin
    dO = -O*d + b*(1 -Di -Dm) + 2*ep*(1 -Di -Dm) - O^2*a*Di
    dDi = -Di * gama + ep*(1 -Di -Dm) + th*Dm + (-1/2)*O^2*a*Di
    dDm = Di * gama - th*Dm
end a ep b gama th d



function first_hit_ss(sol_rd,tspan, eps)
    tss_set = []
    for ti = tspan[1] : 0.05 : tspan[2]
        ϵ = norm(sol_rd(ti) .- sol_rd[end])
        # @show ϵ, ti, sol_rd(ti)
        if ϵ < eps
            push!(tss_set,ti)
             break
        else
            continue
        end
    end
    return tss_set#, sol_rd(ti)
end

# --- Point_T_scaler: Equlibrate time to attractor
p = [1, 0.2, 1., 0.25, 0.25, 0.22]
function Point_T_scaler(x,y,z; tspan = (0.0,1.5e2), p = p)
    u0 = [x,y,z]
    prob_rd = ODEProblem(rn1d_ode,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)

    t_ss = first_hit_ss(sol_rd,tspan,1e-2)[1]
    @show t_ss
    sss = sol_rd[end]
    @show sss
    # attractor_low = norm(sss .- [0,0.05,0]) < norm(sss .- [0.69, 0.73, 0.68])  ? t_ss : 0
    # attractor_high = norm(sss .- 1) > norm(sss .- 0) ? t_ss : 0
    # return attractor_low, attractor_high
    attractor_low  = norm(sss .- [4.3, 0.02, 0.02]) > norm(sss .- [0., 0.5, 0.5]) ? t_ss : 0
    attractor_high = norm(sss .- [4.3, 0.02, 0.02]) < norm(sss .- [0., 0.5, 0.5]) ? t_ss : 0
    return attractor_low, attractor_high
end

Point_T_scaler(0.69,0.7,0.68)
Point_T_scaler(0.3,0.5,0.5)

N = 30 ;r = LinRange(0., 5.5, N)


@time Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])


s1, v1 = textslider(sort(unique(Tc_high))[2]:0.1:sort(unique(Tc_high))[end], "slice")
density1 = Makie.volume(r,r,r,Tc_high, algorithm = :iso, isovalue = v1, color = :green)
scatter!(density1, Point3f0[ (4.3, 0.02, 0.02)], marker = [ :x])
vbox(s1, density1)

s2, v2 = textslider(sort(unique(Tc_low))[2]:0.1:sort(unique(Tc_low))[end], "slice")
density2 = Makie.volume(r,r,r,Tc_low, algorithm = :iso, isovalue = v2, color = :blue, transparency = true)
scatter!(density2, Point3f0[(0., 0.5, 0.5)], marker = [:circle])

vbox(s1, density1, s2, density2)







# Detail look near high_SS (4.3, 0.02, 0.02)

N = 30 ;r = LinRange(0., 6., N)# ;r1 = LinRange(0., 6., N)
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])
s1, v1 = textslider(sort(unique(Tc_high))[2]:0.1:sort(unique(Tc_high))[end], "slice")
density1 = Makie.volume(r,r,r,Tc_high, algorithm = :iso, isovalue = v1, color = :green)
scatter!(density1, Point3f0[ (4.3, 0.02, 0.02),[0.18,0.04,0.48]], marker = [ :circle,:x], markersize = .5f0)
vbox(s1, density1)







# Detail look near low_SS (0.0,0.5,0.5)
N = 30 ;r = LinRange(0., 0.6, N)
@time Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
s2, v2 = textslider(sort(unique(Tc_low))[2]:0.01:sort(unique(Tc_low))[end], "slice")
density2 = Makie.volume(r,r,r,Tc_low, algorithm = :iso, isovalue = v2, color = :blue, transparency = true)
scatter!(density2, Point3f0[(0.,0.5,0.5)], marker = [:circle], markersize = 0.1f0)
vbox(s2, density2)




# ====== volume =====
s, value = textslider(sort(unique(Tc_low))[2]:0.01:sort(unique(Tc_low))[end], "isovalue")
density = Makie.volume(Tc_low, algorithm = :iso, isovalue = value, color = :blue)
hbox(s, density)



# ------------------------------------------------------------------------------------------
# Eliminate Di ------------------

#  ----------- 3D vectorfield ----------
using Parameters, Makie, GeometryTypes,AbstractPlotting
using LinearAlgebra
@with_kw struct rn1d_p{T}
    a::T
    ep::T
    b::T
    gama::T
    th::T
    d::T
end
P = rn1d_p(1., 0.2, 1., 0.25, 0.25, 0.22)
@unpack a, ep, b, gama, th, d = P
# fO  = -O*d + b*Da + 2*ep*Da - O^2*a*(1 -Da -Dm)
# fDa = -ep*Da + (1/2)*O^2*a*(1 -Da -Dm)
# fDm = (1 -Da -Dm)*gama - th*Dm

f(O,Da,Dm,P::rn1d_p) = Point3f0(
    -O*P.d + P.b*Da + 2*P.ep*Da - O^2*P.a*(1 -Da -Dm),
    -P.ep*Da + (1/2)*O^2*P.a*(1 -Da -Dm),
    (1 -Da -Dm)*P.gama - P.th*Dm)
f(x) = f(x, P)


N = 20
r = LinRange(0., 5.5, N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))

vectorfield_dm = f.(x, y, z, Ref(P))
velocity_dm = norm.(vectorfield_dm)
normies_dm = sqrt.(velocity_dm) ./ 10.0
s1, a = textslider(extrema(normies_dm)[1]:0.1:extrema(normies_dm)[2], "Green")
s2, b = textslider(extrema(normies_dm)[1]:0.1:extrema(normies_dm)[2], "Yellow")
s3, c = textslider(extrema(normies_dm)[1]:0.1:extrema(normies_dm)[2], "Red")
vfield_dm = meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield_dm),
    marker = Pyramid(Point3f0(0), .6f0, .6f0),
    markersize = Vec3f0.(0.2, 0.2, vec(normies_dm)),
    color = vec(normies_dm),
    # colormap = lift((x,y,z)-> [(:green, x), (:yellow, y), (:red, z)], a,b,c),
    colormap = [(:green, 0.1), (:yellow, 0.5), (:red, 1.5)],
    colorrange = (extrema(normies_dm)[1], extrema(normies_dm)[2]),
    shading = false,
)
vbox(vfield_dm, colorlegend(vfield_dm[end]))
vbox(s1,s2,s3,vfield_dm, colorlegend(vfield_dm[end]))



# ------------- time to Equilibrate ------------
using DifferentialEquations, ParameterizedFunctions

rn1d_ode = @ode_def begin
    dO = -O*d + b*Da + 2*ep*Da - O^2*a*(1 -Da -Dm)
    dDa = -ep*Da + (1/2)*O^2*a*(1 -Da -Dm)
    dDm = (1 -Da -Dm)*gama - th*Dm
end a ep b gama th d



function first_hit_ss(sol_rd,tspan, eps)
    tss_set = []
    for ti = tspan[1] : 0.05 : tspan[2]
        ϵ = norm(sol_rd(ti) .- sol_rd[end])
        # @show ϵ, ti, sol_rd(ti)
        if ϵ < eps
            push!(tss_set,ti)
             break
        else
            continue
        end
    end
    return tss_set#, sol_rd(ti)
end

# --- Point_T_scaler: Equlibrate time to attractor
p = [1, 0.2, 1., 0.25, 0.25, 0.22]
function Point_T_scaler(x,y,z; tspan = (0.0,1.5e3), p = p)
    u0 = [x,y,z]
    prob_rd = ODEProblem(rn1d_ode,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)

    t_ss = first_hit_ss(sol_rd,tspan,1e-2)[1]
    @show t_ss
    sss = sol_rd[end]
    @show sss

    attractor_low  = norm(sss .- [4.36, 0.96, 0.02]) > norm(sss .- [0., 0., 0.5]) ? t_ss : 0
    attractor_high = norm(sss .- [4.36, 0.96, 0.02]) < norm(sss .- [0., 0., 0.5]) ? t_ss : 0
    return attractor_low, attractor_high
end

Point_T_scaler(0.69,0.7,0.68)
Point_T_scaler(0.,0.,0.5)

N = 50 ;r = LinRange(0., 5.5, N)


@time Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])


s1, v1 = textslider(sort(unique(Tc_high))[2]:0.1:sort(unique(Tc_high))[end], "slice")
density1 = Makie.volume(r,r,r,Tc_high, algorithm = :iso, isovalue = v1, color = :green)
Makie.scatter!(density1, Point3f0[ (4.36, 0.96, 0.02)], marker = [ :x])
vbox(s1, density1)

s2, v2 = textslider(sort(unique(Tc_low))[2]:0.1:sort(unique(Tc_low))[end], "slice")
density2 = Makie.volume(r,r,r,Tc_low, algorithm = :iso, isovalue = v2, color = :blue, transparency = true)
Makie.scatter!(density2, Point3f0[(0., 0., 0.5)], marker = [:circle])

vbox(s1, density1, s2, density2)






# u=unique(Tc_low)
# d=Dict([(i,count(x->x==i,Tc_low)) for i in u])
