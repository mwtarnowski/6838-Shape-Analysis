using LinearAlgebra
using Makie

# Sample Lissajous curve
a = 4
b = 2
delta = pi / 3
n = 100

t = range(0, stop=2*pi, length=n)
x = sin.(a * t * delta)
y = sin.(b * t)
xy = [x y]

## Problem 2.c
n = length(x)
u = zeros(n - 2)
v = zeros(n - 2)
### YOUR CODE HERE TO COMPUTE GRADIENT ###
### END HOMEWORK PROBLEM ###

# Plot curve and gradient vectors
scene1 = Scene(show_axis=false, scale_plot=false)
lines!(scene1, x, y,
       linewidth=3, color=RGBAf0(0.8, 0.05, 0.1))
arrows!(scene1, x[2:end-1], y[2:end-1], u, v,
        arrowsize=0.05)

## Problem 2.d
kappa = zeros(n-2)
### YOUR CODE HERE TO COMPUTE KAPPA ###
### END HOMEWORK PROBLEM ###
c = to_colormap(kappa)

scene2 = Scene(show_axis=false, scale_plot=false)
lines!(scene2, x, y, linewidth=4, color=c)
lines!(scene2,
       zeros(length(c)) .+ (maximum(x) + 0.2),
       range(minimum(y), stop=maximum(y), length=length(c)),
       color=[minimum(c), maximum(c)], linewidth=30)

## Problem 2.e
function topoint2f(xy)
    return Point2f0.(xy[:,1], xy[:,2])
end

t0 = 0
t1 = pi * 1.25
nsamples = 100
# Modify nsteps if your method does not converge
nsteps = 200

# We provide a few examples of curves to try
# curveFunction(t) = [cos(t)-cos(3*t).^3 sin(t)-sin(3*t).^3]
# curveFunction(t) = [cos(t) sin(t)]
curveFunction(t) = [t (t.-t0).*(t1.-t)]
curve = vcat(curveFunction.(range(t0, stop=t1, length=nsamples))...)

scene3 = Scene(show_axis=false, scale_plot=false)
lines!(scene3, curve[:,1], curve[:,2],
       linewidth=3, color=RGBAf0(0.8, 0.05, 0.1))
lineplt = scene3[end]
record(scene3, "straighten.gif", 1:nsteps) do i
    ### YOUR CODE HERE TO PERFORM GRADIENT DESCENT ###
    ### END HOMEWORK PROBLEM ###
    lineplt[1] = curve[:,1]
    lineplt[2] = curve[:,2]
end

# Uncomment to display 2.c/2.d/2.e
# display(scene1)
# display(scene2)
# display(scene3)
