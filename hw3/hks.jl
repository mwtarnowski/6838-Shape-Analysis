using Makie

include("utils.jl")

filename = "data/human_coarse.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)

### ADD CODE TO COMPUTE THE SPECTRUM OF THE LAPLACIAN HERE ###
function laplacianSpectrum(X, T, k)
    return zeros(size(X,1)), rand(size(X,1), k)
end
### END HOMEWORK ASSIGNMENT ###

# Compute 10 eigenfunctions and display eigenfunction number 4
k = 10
vals, vecs = laplacianSpectrum(X, T, k)
scene1 = showdescriptor(X, T, vecs[:,4])

### ADD CODE TO COMPUTE THE HEAT KERNEL SIGNATURE HERE ###
function HKS(eigenvalues, eigenvectors, nSamples)
    neig = size(eigenvalues, 1)
    tmin = 4*log(10)/eigenvalues[neig]
    tmax = 4*log(10)/eigenvalues[2]
    t = exp10.(range(tmin, stop=tmax, length=nSamples))

    hks = rand(size(eigenvectors, 1), nSamples)

    return hks
end
### END HOMEWORK ASSIGNMENT ###

nSamples = 500
k = 300

vals, vecs = laplacianSpectrum(X, T, k)
hks = HKS(vals, vecs, nSamples)

hand1, hand2 = 259, 135
knee1, knee2 = 232, 257

x = range(1, stop=nSamples, length=nSamples)

fig = Figure()
trace0 = scatter(fig[1,1], x, hks[hand1,:], markersize=5)
trace1 = scatter(fig[1,2], x, hks[hand2,:], markersize=5)
trace2 = scatter(fig[2,1], x, hks[knee1,:], markersize=5)
trace3 = scatter(fig[2,2], x, hks[knee2,:], markersize=5)
fig

### COMPUTE DIFFERENCE FUNCTION |HKS(x0) - HKS(x)|_2 HERE ###
result = rand(nv)
### END HOMEWORK ASSIGNMENT ###
scene3 = showdescriptor(X, T, result)
