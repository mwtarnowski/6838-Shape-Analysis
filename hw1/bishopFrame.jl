using Makie
using LinearAlgebra

function topoint3f(xyz)
    return Point3f0.(xyz[:,1], xyz[:,2], xyz[:,3])
end

# Sample Lissajous curve
a = 3
b = 2
delta = pi / 2
n = 100

t = range(0, stop=2*pi, length=n)
x = sin.(a * t * delta)
y = sin.(b * t)
z = .5 * t
xyz = [x y z]

# Edges and midpoints
diff = xyz[2:end,:] - xyz[1:end-1,:]
midpoints = (xyz[1:end-1,:] + xyz[2:end,:]) / 2

## Problem 3.a
binormal = zeros(n-2,3)
### YOUR CODE TO COMPUTE BINORMAL HERE
### END HOMEWORK PROBLEM ###

scene1 = Scene()
lines!(scene1, x, y, z,
       linewidth=4, color=:black)
arrows!(scene1, topoint3f(xyz[2:end-1,:]), topoint3f(binormal),
        linecolor=:red, arrowcolor=:red,
        linewidth=2, arrowsize=0.03)

## Problem 3.b
### WRITE ANY PRE-COMPUTATION CODE HERE ###
### END PRE-COMPUTATION ###

scene2 = Scene()
lines!(scene2, x, y, z,
       linewidth=4, color=:black)
arrows!(scene2, topoint3f(midpoints),
        topoint3f(zeros(n-1, 3)),
        linecolor=:red, arrowcolor=:red,
        linewidth=2, arrowsize=0.03)
arrows!(scene2, topoint3f(midpoints),
        topoint3f(zeros(n-1, 3)),
        linecolor=:blue, arrowcolor=:blue,
        linewidth=2, arrowsize=0.03)
uplot = scene2[end-1]
vplot = scene2[end]

scale = 0.25
u = zeros(n - 1, 3)
v = zeros(n - 1, 3)
record(scene2, "bishop.gif", 0:5:360) do angle
    t = diff[1,:] / norm(diff[1,:])
    zaxis = [0;0;1]
    frame0 = cross(t, zaxis)
    frame1 = cross(frame0, t)

    theta = angle * pi / 180
    u[1,:] = cos(theta) .* frame0 + sin(theta) .* frame1
    v[1,:] = -sin(theta) .* frame0 + cos(theta) .* frame1
    for frame in 2:(n-1)
        u[frame,:] .= 0
        v[frame,:] .= 0
        ### YOUR CODE TO FILL IN u AND v HERE ###
        ### END HOMEWORK PROBLEM ###
    end
    uplot[2] = topoint3f(u.*scale)
    vplot[2] = topoint3f(v.*scale)
end


## Problem 3.c
# thetas is the variable you are supposed to update.
thetas = range(0, stop=3*pi, length=n-1)
m1 = cos.(thetas).*u + sin.(thetas).*v
m2 = -sin.(thetas).*u + cos.(thetas).*v

scene3 = Scene()
lines!(scene3, x, y, z,
       linewidth=4, color=:black)
arrows!(scene3, topoint3f(midpoints),
        topoint3f(scale.*m1),
        linecolor=:red, arrowcolor=:red,
        linewidth=2, arrowsize=0.03)
arrows!(scene3, topoint3f(midpoints),
        topoint3f(scale.*m2),
        linecolor=:blue, arrowcolor=:blue,
        linewidth=2, arrowsize=0.03)
m1plot = scene3[end-1]
m2plot = scene3[end]

nsteps = 200
grad = zeros(n-1)
record(scene3, "untwist.gif", 1:nsteps) do i
    ### YOUR CODE HERE TO UPDATE thetas VIA GRADIENT DESCENT ###
    ### END HOMEWORK PROBLEM ###
    m1 = cos.(thetas).*u + sin.(thetas).*v
    m2 = -sin.(thetas).*u + cos.(thetas).*v
    m1plot[2] = topoint3f(scale.*m1)
    m2plot[2] = topoint3f(scale.*m2)
end

# Uncomment to plot scene1
# display(scene1)
