using Makie


# =========== 1. 3D streamplot example
struct FitzhughNagumo{T}
    ϵ::T; s::T; γ::T; β::T
end

P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)
f(x, P::FitzhughNagumo) = Point3f0(
    (x[1]-x[2]-x[1]^3+P.s)/P.ϵ,
    P.γ*x[2]-x[2] + P.β,
    P.γ*x[1]-x[3] - P.β,)
f(x) = f(x, P) # why need this step? 🔯

f(position, P::FitzhughNagumo) = Point3f0(
    let x, y, z = position,
    (x-y-x^3+P.s)/P.ϵ,
    P.γ*y-y + P.β,
    P.γ*x-z - P.β,)
streamplot(f, -1.5..1.5, -1.5..1.5, -1.5..1.5, colormap = :magma, gridsize = (10, 10), arrow_size = 0.06)

# scene = Scene()
# s = streamplot(scene,f, -1.5..1.5, -1.5..1.5, -1.5..1.5, colormap = :magma, gridsize = (10, 10), arrow_size = 0.06)







# =========== 2. GUI for exploring Lorenz equation
using Makie
using Colors
using AbstractPlotting: textslider, colorswatch

s1, a = textslider(0f0:50f0, "a", start = 13)
s2, b = textslider(-20f0:20f0, "b", start = 10)
s3, c = textslider(0f0:20f0, "c", start = 2)
s4, d = textslider(range(0.0, stop = 0.02, length = 100), "d", start = 0.01)
s5, scales = textslider(range(0.01, stop = 0.5, length = 100), "scale", start = 0.1)
s6, colorsw, pop = colorswatch()
# Lorenz function
function lorenz(t0, a, b, c, h)
 Point3f0(
     t0[1] + h * a * (t0[2] - t0[1]),
     t0[2] + h * (t0[1] * (b - t0[3]) - t0[2]),
     t0[3] + h * (t0[1] * t0[2] - c * t0[3]),
 )
end
# step through the `time`
function lorenz(array::Vector, a = 5.0 ,b = 2.0, c = 6.0, d = 0.01)
 t0 = Point3f0(0.1, 0, 0)
 for i = eachindex(array)
     t0 = lorenz(t0, a,b,c,d)
     array[i] = t0
 end
 array
end

n1, n2 = 18, 30
N = n1*n2
args_n = (a, b, c, d)
v0 = lorenz(zeros(Point3f0, N), to_value.(args_n)...)
positions = lift(lorenz, Node(v0), args_n...)
rotations = lift(diff, positions)
rotations = lift(x-> push!(x, x[end]), rotations)
mesh_scene = meshscatter(
positions,
markersize = scales, rotation = rotations,
intensity = collect(range(0f0, stop = 1f0, length = length(positions[]))),
color = colorsw
)
# parent = Scene(resolution = (1000, 800))
vbox(
hbox(s1, s2, s3, s4, s5, s6),
mesh_scene
)







# ============ 3. Simple contour plot
# contour plot
function f(x,y,z)
    v1 = x-y-x^3
    v2 = y^2-y + 29
    v3 = x-z - y*z
    return v1^2 + v2^2 + v3^2
end

rg = LinRange(-1.5, 1.5, 10)
 r = range(-2, stop = 2, length = 100)
me = [f(x,y,z) for x = r ,y =r, z=r]
contour(me)







# ======= 2d streamplot =========
using Makie
using Random
# struct FitzhughNagumo{T}
#     ϵ::T
#     s::T
#     γ::T
#     β::T
# end
Random.seed!(10)
point(x) = Makie.Point2f0(x)
point(x, y) = Makie.Point2f0(x, y)
const Point = Makie.Point2f0

# f(x, P::FitzhughNagumo) = Point((x[1]-x[2]-x[1]^3+P.s)/P.ϵ,
#                                          P.γ*x[1]-x[2] + P.β)
#
# f(x) = f(x, P)

# new f(x) with DC model
f(x) = Point( (0.276*x[1]^2 + 1.38*x[2]^2 + 0.897*x[1]^2*x[2]^2)/(1+ x[1]^2 + x[2]^2 + x[1]^2*x[2]^2) - 0.138*x[1],
(0.00828*x[1]^2 + 0.0828*x[2]^2 + 0.092*x[1]^2*x[2]^2)/(1+ x[1]^2 + x[2]^2 + x[1]^2*x[2]^2) - 0.046*x[2] )

# P = FitzhughNagumo(0.1, 0.0, 1.5, 0.8)
R = 1.5
l = (32, 32)
limits = FRect(0, 0, 10R, 10R)
r = (range(limits.origin[1], limits.origin[1] + limits.widths[1], length=l[1]+1),
 range(limits.origin[2], limits.origin[2] + limits.widths[2], length=l[2]+1)
)
p = Scene(resolution=(1600,1000), limits=limits)

mask = trues(l...)
arrows = Tuple{Point,Point}[]

for c in CartesianIndices(mask)

    i0, j0 = Tuple(c)
    x0 = Point(first(r[1]) + (i0 - 0.01)*step(r[1]),
              first(r[2]) + (j0 - 0.01)*step(r[1]))
    dt = 0.005

    if mask[c]
        push!(arrows, (x0, f(x0)/norm(f(x0))))

        mask[c] == false
        for d in (-1,1)
            x = x0
            xx = [x]
            i0, j0 = Tuple(c)
            while x in limits && length(xx) < 500
                x = x + d*dt*f(x)/norm(f(x))
                if !(x in limits)
                    break
                end
                i = searchsortedlast(r[1], x[1])
                j = searchsortedlast(r[2], x[2])
                if (i,j) != (i0, j0)
                    if !mask[i,j]
                        break
                    end
                    mask[i,j] = false
                    i0, j0 = i, j
                end
                push!(xx, x)
            end
            lines!(p, xx)
        end
    end

end
arrows!(p, first.(arrows), last.(arrows), lengthscale=0.0, arrowsize=0.3)
display(p)
