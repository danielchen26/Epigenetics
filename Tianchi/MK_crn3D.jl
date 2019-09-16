using Makie
using AbstractPlotting
using GeometryTypes,LinearAlgebra


struct FitzhughNagumo{T}
    ϵ::T; s::T; γ::T; β::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)
# f(x, P::FitzhughNagumo) = Point3f0(
#     (x[1]-x[2]-x[1]^3+P.s)/P.ϵ,
#     P.γ*x[2]-x[2] + P.β,
#     P.γ*x[1]-x[3] - P.β,)

f(x,y,z, P::FitzhughNagumo) = Point3f0(
    (x-y-x^3+P.s)/P.ϵ,
    P.γ*y-y + P.β,
    P.γ*x-z - P.β,)
# f(x) = f(x, P)

# P = lift(FitzhughNagumo,to_value(args_n)...)
# f(x) = f(x,P)
# s_plot = streamplot(f, -1.5..1.5, -1.5..1.5, -1.5..1.5, colormap = :magma, gridsize = (10, 10), arrow_size = 0.06)
# hbox(s,s_plot)#


N = 20
r = LinRange(-1.5, 1.5, N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))

vectorfield_dm = f.(x, y, z, Ref(P))
velocity_dm = norm.(vectorfield_dm)
normies_dm = sqrt.(velocity_dm) ./ 10.0
s1, a = textslider(0.07:0.1:0.58, "Green")
s2, b = textslider(0.07:0.1:0.58, "Yellow")
s3, c = textslider(0.07:0.1:0.58, "Red")
vfield_dm = meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield_dm),
    marker = Pyramid(Point3f0(0), 1f0, 1f0),
    markersize = Vec3f0.(0.02, 0.02, vec(normies_dm)),
    color = vec(normies_dm),
    # colormap = lift((x,y,z)-> [(:green, x), (:yellow, y), (:red, z)], a,b,c),
    colormap = [(:green, extrema(normies_dm)[1]), (:yellow, extrema(normies_dm)[2]), (:red, extrema(normies_dm)[2]+1)],
    colorrange = (extrema(normies_dm)[1], extrema(normies_dm)[2]),
    shading = false,
)
vbox(vfield_dm, colorlegend(vfield_dm[end]))
vbox(s1,s2,s3,vfield_dm, colorlegend(vfield_dm[end]))











# ------------------------------------------------

using Makie
using LinearAlgebra

struct FitzhughNagumo{T}
    ϵ::T; s::T; γ::T; β::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)

f(x,y,z, P::FitzhughNagumo) = Point3f0(
    (x-y-x^3+P.s)/P.ϵ,
    P.γ*y-y + P.β,
    P.γ*x-z - P.β,)

N = 20
r = LinRange(-1.5, 1.5, N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))

vectorfield_dm = f.(x, y, z, Ref(P))
velocity_dm = norm.(vectorfield_dm)
normies_dm = sqrt.(velocity_dm) ./ 10.0
s1, a = textslider(0:0.1:0.6, "Green")
s2, b = textslider(0:0.1:0.6, "Yellow")
s3, c = textslider(0:0.1:0.6, "Red")
vfield_dm = meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield_dm),
    marker = Pyramid(Point3f0(0), 1f0, 1f0),
    markersize = Vec3f0.(0.02, 0.02, vec(normies_dm)),
    color = vec(normies_dm),
    colormap = lift((x,y,z)-> [(:green, x), (:yellow, y), (:red, z)], a,b,c),
    shading = false,
)
hbox(s1,s2,s3,vfield_dm)
# ------------------------------------------------

vbox(s1, s2, s3, vfield_dm, colorlegend(vfield_dm[end]))































#  ------------- My Reduced model example -------------
using Parameters, Makie
using AbstractPlotting
@with_kw struct rd{T}
    Kₒ::T
    Kₙₜ::T
    Kd::T
    a₁::T
    d::T
    aₙₜ::T
    aₒ::T
    αₜ::T
    αₒ::T
    αₙ::T
    δ::T
    β::T
    m1::T
    m2::T
    m3::T
end
P = rd(0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6)
@unpack Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3 = P

# =====  3D streamplot ========
f(N,T,O) = Point3f0(
    -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1,

    (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))),

    (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    )
streamplot(f, -1.5..1.5, -1.5..1.5, -1.5..1.5, colormap = :magma, gridsize = (10, 10), arrow_size = 0.04)

# ------ volume control
using Parameters, Makie
@with_kw struct rd{T}
    Kₒ::T
    Kₙₜ::T
    Kd::T
    a₁::T
    d::T
    aₙₜ::T
    aₒ::T
    αₜ::T
    αₒ::T
    αₙ::T
    δ::T
    β::T
    m1::T
    m2::T
    m3::T
end
P = rd(0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6)
@unpack Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3 = P
function f(N,T,O)
    v1 = -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1

    v2 = (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))

    v3 = (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    return sqrt(v1^2 + v2^2 + v3^2)
end
r = range(0, stop = 1, length = 100)
me = [f(x,y,z) for x = r ,y =r, z=r]
Makie.contour(r,r,r, log.(me), levels = 6, alpha = 0.3, transparency = true)

# ------ heatmap control
s, value = textslider(1:size(me, 3), "slice")
hmap = heatmap(lift(idx-> me[:, :, idx], value), padding = (0.0, 0.0))
hbox(s, hmap)

# ------ 2D contour control
s, value = textslider(1:size(me, 3), "slice")
contour_2d = Makie.contour(lift(idx-> me[:, :, idx], value), levels = 10, colormap = :viridis, padding = (0.0, 0.0))
hbox(s, contour_2d)

velo_log = log.(me)
s, value = textslider(LinRange(extrema(velo_log)..., 100), "isovalue")
density = Makie.volume(velo_log, algorithm = :iso, isovalue = value, color = :blue)
hbox(s, density)

density2 = contour(velo_log, algorithm = :iso, isovalue = value,levels = reverse([2.1, -0.7, -2.5, -2.97, -3.7]), alpha = 0.2)
hbox(s, density2)

contour(r, r,r, velo_log, levels = reverse([2.1, -0.7, -2.5, -2.97, -3.7]), alpha = 0.2)




# add textslider to streamplot

using Parameters, Makie
@with_kw struct rd{T}
    Kₒ::T
    Kₙₜ::T
    Kd::T
    a₁::T
    d::T
    aₙₜ::T
    aₒ::T
    αₜ::T
    αₒ::T
    αₙ::T
    δ::T
    β::T
    m1::T
    m2::T
    m3::T
end
P = rd(0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6)
@unpack Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3 = P
function f(N,T,O)
    v1 = -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1

    v2 = (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))

    v3 = (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    return sqrt(v1^2 + v2^2 + v3^2)
end
r = range(-2, stop = 2, length = 100)
me = [f(x,y,z) for x = r ,y =r, z=r]
contour(r,r,r, log.(me), levels = 6, alpha = 0.3, transparency = true)

s, Kd_v = textslider(1:size(Kd, 3), "slice")














# ======== to play with vectorfield ========
using Parameters, Makie
using AbstractPlotting
using GeometryTypes,LinearAlgebra

@with_kw struct rd{T}
    Kₒ::T
    Kₙₜ::T
    Kd::T
    a₁::T
    d::T
    aₙₜ::T
    aₒ::T
    αₜ::T
    αₒ::T
    αₙ::T
    δ::T
    β::T
    m1::T
    m2::T
    m3::T
end
P = rd(0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6)
@unpack Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3 = P

f(N,T,O) = Point3f0(
    -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1,

    (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))),

    (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    )

# x = range(-2, stop = 2, length = 100)

N = 50
r = LinRange(0, 1., N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))
vectorfield = f.( x, y, z)
velocity = norm.(vectorfield)
normies = sqrt.(velocity) ./ 10.0


s1, a = textslider(0:0.1:1., "Green")
s2, b = textslider(0:0.1:1., "Yellow")
s3, c = textslider(0:0.1:1., "Red")

vfield = meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield),
    marker = Pyramid(Point3f0(0), 0.5f0, 0.5f0),
    markersize = Vec3f0.(0.04, 0.04, vec(normies)),
    color = vec(normies),
    colormap = lift((cv1,cv2,cv3)-> [(:green, cv1), (:yellow, cv2), (:red, cv3)], a,b,c),
    shading = false,
)
hbox(s1,s2,s3,vfield)













# ======== to play with vectorfield another test following above========
using Parameters, Makie

@with_kw struct rd{T}
    Kₒ::T
    Kₙₜ::T
    Kd::T
    a₁::T
    d::T
    aₙₜ::T
    aₒ::T
    αₜ::T
    αₒ::T
    αₙ::T
    δ::T
    β::T
    m1::T
    m2::T
    m3::T
end
P = rd(0.3, 0.2, 0.1, 1., 1., 1000., 1000., 1.0, 1.0, 1.0, 1., 1e-6, 1e-6, 0.05, 1e-6)
@unpack Kₒ, Kₙₜ, Kd, a₁, d, aₙₜ, aₒ, αₜ, αₒ, αₙ, δ, β, m1, m2, m3 = P
f(N,T,O) = Point3f0(
    -Kd*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - N*T*a₁ - N*δ + O*αₙ/(Kₒ + O) + m1,

    (-Kd*Kₒ*Kₙₜ*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*Kₒ*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kd*Kₙₜ*O*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)/(2*d*β) - Kd*O*a₁*(d + β)^2*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^3/(16*d^2*β^3) - Kₒ*Kₙₜ*N*T*a₁ - Kₒ*Kₙₜ*T*δ + Kₒ*Kₙₜ*m2 - Kₒ*N*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₒ*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*N*O*T*a₁ - Kₙₜ*O*T*δ + Kₙₜ*O*m2 - N*O*T*a₁*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - O*T*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m2*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₜ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))),

    (-Kₒ*Kₙₜ*O*δ + Kₒ*Kₙₜ*m3 - Kₒ*O*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + Kₒ*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) - Kₙₜ*O^2*δ + Kₙₜ*O*m3 - O^2*δ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*m3*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2) + O*αₒ*(d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2))/((Kₒ + O)*(Kₙₜ + (d + β)*(Kd*a₁ + β - ((Kd^2*a₁^2*d + Kd^2*a₁^2*β + 2*Kd*a₁*d*β + 2*Kd*a₁*β^2 + 4*N*T*a₁*d*β + d*β^2 + β^3)/(d + β))^0.5)^2/(8*d*β^2)))
    )

using AbstractPlotting
using GeometryTypes,LinearAlgebra
# x = range(-2, stop = 2, length = 100)

N = 20
r = LinRange(0, 1., N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))

vectorfield = f.( x, y, z)
velocity = norm.(vectorfield)
# normies = sqrt.(velocity) ./ 10.0

s1, v1 = textslider(0:0.05:1.5, "slice")
s2, v2 = textslider(0:0.05:1.5, "slice")
s3, v3 = textslider(0:0.05:1.5, "slice")
# value = sort(unique(velocity))[1:1000:end]






# --- manually change the norm velocity -------
v1 = 0.2
n_velocity = velocity .*(to_value(v1) .< velocity .< to_value(v1) .+ 0.1)

# f(a) = velocity .*(a .< velocity .< a .+ 0.1)
# n_velocity = lift(a -> f(a), 0.1)
n_velocity =Tc_high./10
vfield = Makie.meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield),
    marker = Pyramid(Point3f0(0), 1f0, 1f0), # middle::Point{3, T}, length::T, width ::T
    markersize = Vec3f0.(0.2, 0.2, vec(n_velocity./10)),
    color = vec(n_velocity .+0.01),
    colormap = :Spectral, # [(:green, 0.1), (:yellow, 0.5), (:red, 1.5)],
    shading = false,
)
hbox(s1,vfield)




# ---- slider control the contour with iso value that define the velocity ----- ready to show
s1, v1 = textslider(0.001:0.01:1.5, "slice")
density = Makie.volume(r,r,r,velocity, algorithm = :iso, isovalue = v1, color = :green)
hbox(s1, density)













# -------- paper_reduce_ode plotting 3d time contour ---------.
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

# ----- for a 3D unit cubit, show SSS --------
# N = 1; sp = 0.1; r = 0:sp:N
# tspan = (0.,200.)
# scene = Scene(resolution = (500, 500))
# for x = r ,y =r, z=r
#     u0 = [x,y,z]
#     println("Initial :", u0)
#     prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
#     sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
#     println("SSS :\n", sol_rd[end])
#     println("Equilibrate Time:", first_hit_ss(sol_rd,tspan, 1e-2)[1])
# end
N = 30 ;r = LinRange(0., 1.5, N)
# x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))
# positions = vec(Point3f0.(x, y, z))
# scene = Scene(resolution = (500, 500))
# v1 = 1
# t_MF = Tc_high .*(to_value(v1) .< Tc_high .< to_value(v1) .+ 2)


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
function Point_T_scaler(x,y,z; tspan = (0.0,1.5e2), p = p)
    u0 = [x,y,z]
    prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)

    t_ss = first_hit_ss(sol_rd,tspan,1e-3)[1]
    @show t_ss
    sss = sol_rd[end]
    # @show sss
    # attractor_low = norm(sss .- [0,0.05,0]) < norm(sss .- [0.69, 0.73, 0.68])  ? t_ss : 0
    # attractor_high = norm(sss .- 1) > norm(sss .- 0) ? t_ss : 0
    # return attractor_low, attractor_high
    attractor_low  = norm(sss .- 0) < norm(sss .- [0.69, 0.73, 0.68]) ? t_ss : 0
    attractor_high = norm(sss .- [0.69, 0.73, 0.68]) < norm(sss .- 0) ? t_ss : 0
    return attractor_low, attractor_high
end

Point_T_scaler(0.69,0.7,0.68)
Point_T_scaler(0.,0.03,0.0)

@time Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])


s1, v1 = textslider(sort(unique(Tc_high))[2]:0.1:sort(unique(Tc_high))[end], "slice")
density1 = Makie.volume(r,r,r,Tc_high, algorithm = :iso, isovalue = v1, color = :green)
scatter!(density1, Point3f0[ (0.695,0.734,0.685), [0.3,0.18,0.13]], marker = [ :circle, :x])
# vbox(s1, density1)

s2, v2 = textslider(sort(unique(Tc_low))[2]:0.1:sort(unique(Tc_low))[end], "slice")
density2 = Makie.volume(r,r,r,Tc_low, algorithm = :iso, isovalue = v2, color = :blue, transparency = true)
scatter!(density2, Point3f0[(0,0.05,0), [0.3,0.18,0.13]], marker = [:circle,:x])

vbox(s1, density1, s2, density2)




# Detail look near high_SS (0.695,0.734,0.685)

N = 30 ;r = LinRange(0.6, 0.8, N)
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
s1, v1 = textslider(sort(unique(Tc_high))[2]:0.1:sort(unique(Tc_high))[end], "slice")
density1 = Makie.volume(r,r,r,Tc_high, algorithm = :iso, isovalue = v1, color = :green)
scatter!(density1, Point3f0[ (0.695,0.734,0.685)], marker = [ :circle], markersize = 0.01f0)
vbox(s1, density1)







# Detail look near low_SS (0,0.05,0)
N = 30 ;r = LinRange(0., 0.1, N)
@time Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])

s2, v2 = textslider(sort(unique(Tc_low))[2]:0.01:sort(unique(Tc_low))[end], "slice")
density2 = Makie.volume(r,r,r,Tc_low, algorithm = :iso, isovalue = v2, color = :blue, transparency = true)
scatter!(density2, Point3f0[(0,0.05,0)], marker = [:circle], markersize = 0.01f0)
vbox(s2, density2)



# tspan=(0.,200.)
# u0=[0.6,0.7,0.6] # 12
# u0=[0.69,0.73,0.68] # 9
# u0=[0.695,0.734,0.684] # 7
# u0=[0.6953,0.7348,0.6849]
# u0=[0.3,0.3,0.3]
# tspan = (0.,200.)
# prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
# sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
# sol_rd
#
# first_hit_ss(sol_rd,tspan, 1e-2)
#
# using Plots
# plot(sol_rd)
# u0=[0.3,0.3,0.3]
# prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
# sol_rd2 = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
# plot!(sol_rd2)
#
# plotly()








# =========================== Below  is demethylation  ================================
# =========================== Below  is demethylation  ================================
using DifferentialEquations, ParameterizedFunctions, LinearAlgebra, Measures, Queryverse
using DataFrames, Makie, Images #Plots, VegaLite
using Distributions, ProgressMeter, LinearAlgebra
module epi
    include("functions.jl")
    include("./Model_collection/ODE_model_set.jl")
end

# For Demethy_reduced_ode 3d plot

Demethy_reduced_ode_MatlabC = @ode_def_bare begin # m1 -> N     m2 -> T     m3 -> O

    dN = (N*alpha1*r1 + O*alpha1*xi1 + K0*m1*r1 + N*m1*r1 + K0*m1*xi1 + O*m1*xi1 - N^2*delta0*r1 - K0*N*delta0*r1 - K0*N*delta0*xi1 - N*O*delta0*xi1 - K0*N*r1*xi1 + K0*O*r1*xi1)/(K0*r1 + N*r1 + K0*xi1 + O*xi1)

    dT = m2 - T*delta0 + (O*alpha1)/(K0 + O)

    dO = (Kd*(K0*a0*beta*m3 - K0*O*a0*beta*delta0 - K0*O*a0*beta*xi1 + (K0^2*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (K0*O*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (N*O*T*alpha1*beta*k3)/Kd + (K0*N*T*a0*k3*m3)/Kd + (K0*N*T*beta*k3*m3)/Kd + (N*O*T*beta*k3*m3)/Kd - (N*O^2*T*beta*delta0*k3)/Kd - (N*O^2*T*beta*k3*xi1)/Kd - (K0*N*O*T*a0*delta0*k3)/Kd - (K0*N*O*T*beta*delta0*k3)/Kd - (K0*N*O*T*a0*k3*xi1)/Kd - (K0*N*O*T*beta*k3*xi1)/Kd + (K0^2*N*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0^2*N*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (N*O^2*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0*N*O*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (2*K0*N*O*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1))))/(K0*Kd*a0*beta + K0*N*T*a0*k3 + K0*N*T*beta*k3 + N*O*T*beta*k3)
end K0 Kd r1 alpha1 delta0 a0 beta xi1 k3 m1 m2 m3

u0 = rand(3); tspan = (0.,1000.)
p = [0.3, 0.1, 1., 1., 1., 0.5, 1., 1., 1., 0.,0.05,0.]

# ----- for a 3D unit cubit, show SSS --------
# N = 1; sp = 0.1; r = 0:sp:N
# tspan = (0.,200.)
# scene = Scene(resolution = (500, 500))
# for x = r ,y =r, z=r
#     u0 = [x,y,z]
#     println("Initial :", u0)
#     prob_rd = ODEProblem(paper_reduce_ode,u0,tspan,p)
#     sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)
#     println("SSS :\n", sol_rd[end])
#     println("Equilibrate Time:", first_hit_ss(sol_rd,tspan, 1e-2)[1])
# end
N = 30 ;r = LinRange(0., 1.5, N)
# x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))
# positions = vec(Point3f0.(x, y, z))
# scene = Scene(resolution = (500, 500))
# v1 = 1
# t_MF = Tc_high .*(to_value(v1) .< Tc_high .< to_value(v1) .+ 2)

# --- Point_T_scaler: Equlibrate time to attractor
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


function Point_T_scaler(x,y,z; tspan = (0.0,1.5e2), p = p)
    u0 = [x,y,z]
    prob_rd = ODEProblem(Demethy_reduced_ode_MatlabC,u0,tspan,p)
    sol_rd = solve(prob_rd, Rosenbrock23(),abstol=1e-8)

    t_ss = first_hit_ss(sol_rd,tspan,1e-2)[1]
    @show t_ss
    sss = sol_rd[end]
    # @show sss

    attractor_low  = norm(sss .- 0) < norm(sss .- [0.648,0.645,0.546]) ? t_ss : 0
    attractor_high = norm(sss .- [0.648,0.645,0.546]) < norm(sss .- 0) ? t_ss : 0
    return attractor_low, attractor_high
end

Point_T_scaler(0.69,0.7,0.68)
Point_T_scaler(0.,0.2,0.)

@time Tc_low = Array{Float64}([Point_T_scaler(x,y,z)[1] for x = r ,y =r, z=r])
@time Tc_high = Array{Float64}([Point_T_scaler(x,y,z)[2] for x = r ,y =r, z=r])


s1, v1 = textslider(sort(unique(Tc_high))[1]:0.1:sort(unique(Tc_high))[end], "slice")
density1 = Makie.volume(r,r,r,Tc_high, algorithm = :iso, isovalue = v1, color = :green)
scatter!(density1, Point3f0[(0.649, 0.696, 0.548)], marker = [ :circle])
# vbox(s1, density1)

s2, v2 = textslider(sort(unique(Tc_low))[1]:0.1:sort(unique(Tc_low))[end], "slice")
density2 = Makie.volume(r,r,r,Tc_low, algorithm = :iso, isovalue = v2, color = :blue, transparency = true)
scatter!(density2, Point3f0[(0,0.05,0)], marker = [:x])

vbox(s1, density1, s2, density2)









# ======== Demethy_crn vectorfield ========

using Parameters, Makie
@with_kw struct demethy_rd{T}
    K0::T
    Kd::T
    r1::T
    alpha1::T
    delta0::T
    a0::T
    beta::T
    xi1::T
    k3::T
    m1::T
    m2::T
    m3::T
end
P = demethy_rd(0.3, 0.1, 1., 1., 1., 0.5, 1., 1., 1., 0.,0.05,0.)
@unpack K0, Kd, r1, alpha1, delta0, a0, beta, xi1, k3, m1, m2, m3 = P

f_demethy(N,T,O) = Point3f0(
    (N*alpha1*r1 + O*alpha1*xi1 + K0*m1*r1 + N*m1*r1 + K0*m1*xi1 + O*m1*xi1 - N^2*delta0*r1 - K0*N*delta0*r1 - K0*N*delta0*xi1 - N*O*delta0*xi1 - K0*N*r1*xi1 + K0*O*r1*xi1)/(K0*r1 + N*r1 + K0*xi1 + O*xi1),

    m2 - T*delta0 + (O*alpha1)/(K0 + O),

    (Kd*(K0*a0*beta*m3 - K0*O*a0*beta*delta0 - K0*O*a0*beta*xi1 + (K0^2*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (K0*O*a0*beta*xi1*(N*r1 + O*xi1))/(K0*r1 + N*r1 + K0*xi1 + O*xi1) + (N*O*T*alpha1*beta*k3)/Kd + (K0*N*T*a0*k3*m3)/Kd + (K0*N*T*beta*k3*m3)/Kd + (N*O*T*beta*k3*m3)/Kd - (N*O^2*T*beta*delta0*k3)/Kd - (N*O^2*T*beta*k3*xi1)/Kd - (K0*N*O*T*a0*delta0*k3)/Kd - (K0*N*O*T*beta*delta0*k3)/Kd - (K0*N*O*T*a0*k3*xi1)/Kd - (K0*N*O*T*beta*k3*xi1)/Kd + (K0^2*N*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0^2*N*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (N*O^2*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (K0*N*O*T*a0*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1)) + (2*K0*N*O*T*beta*k3*xi1*(N*r1 + O*xi1))/(Kd*(K0*r1 + N*r1 + K0*xi1 + O*xi1))))/(K0*Kd*a0*beta + K0*N*T*a0*k3 + K0*N*T*beta*k3 + N*O*T*beta*k3)
    )

using AbstractPlotting
using GeometryTypes,LinearAlgebra
# x = range(-2, stop = 2, length = 100)

N = 20
r = LinRange(0, 1., N)
x, y, z = reshape(r, (N, 1, 1)), reshape(r, (1, N, 1)), reshape(r, (1, 1, N))

vectorfield_dm = f_demethy.( x, y, z)
velocity_dm = norm.(vectorfield_dm)
normies_dm = sqrt.(velocity_dm) ./ 10.0
vfield_dm = meshscatter(
    vec(Point3f0.(x, y, z)),
    rotation = vec(vectorfield_dm),
    marker = Pyramid(Point3f0(0), 1f0, 1f0),
    markersize = Vec3f0.(0.02, 0.02, vec(normies_dm)),
    color = vec(normies_dm),
    colormap = [(:green, 0.1), (:yellow, 0.5), (:red, 1.0)],
    shading = false,
)




vbox(vfield, vfield_dm)























# ------- Local used function ---- run this first ------

# === The first time solution hit SSS within ϵ tolerance
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












# # Gradient Symbol/Strings:
#     Accent
#     Blues
#     BrBG
#     BuGn
#     BuPu
#     Dark2
#     GnBu
#     Greens
#     Greys
#     OrRd
#     Oranges
#     PRGn
#     Pastel2
#     PiYG
#     PuBu
#     PuBuGn
#     PuOr
#     PuRd
#     Purples
#     RdBu
#     RdGy
#     RdPu
#     RdYlBu
#     RdYlGn
#     Reds
#     Set2
#     Spectral
#     YlGn
#     YlGnBu
#     YlOrBr
#     YlOrRd
#     algae
#     amp
#     balance
#     bgy
#     bgyw
#     bjy
#     bkr
#     bky
#     blues
#     bluesreds
#     bmw
#     cinferno
#     colorwheel
#     coolwarm
#     curl
#     cyclic1
#     cyclic2
#     cyclic3
#     darkrainbow
#     darktest
#     deep
#     delta
#     dense
#     diff
#     dimgreys
#     fire
#     grays
#     grayscale
#     greens
#     greys
#     gwv
#     haline
#     heat
#     ice
#     inferno
#     isolum
#     kb
#     kdc
#     kg
#     kgy
#     kr
#     lightrainbow
#     lighttest
#     magma
#     matter
#     oxy
#     phase
#     plasma
#     pu_or
#     rain
#     rainbow
#     reds
#     redsblues
#     solar
#     speed
#     tarn
#     tempo
#     thermal
#     topo
#     turbid
#     viridis
