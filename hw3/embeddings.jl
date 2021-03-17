using DataFrames
using CSV
using Makie
using Arpack
using NearestNeighbors

include("utils.jl")

filename = "data/swissroll.txt"
data = CSV.read(filename, DataFrame; datarow=1, delim=' ', ignorerepeated=true)
X = [data.Column1 data.Column2 data.Column3]

scene = scatter(X[:,1], X[:,2], X[:,3], color=X[:,3], markersize=0.3)

#### ADD CODE FOR YOUR TWO METHODS HERE ####
n = size(X, 1)

function MDS(D)
    D2 = D .* D
    C = I - ones(n,n) ./ n
    B = -C * D2 * C ./ 2
    vals, vecs = eigs(B, nev=2)
    Y = vecs * diagm(0 => sqrt.(vals))
    return Y
end

function pathDistances(X, K)
    kdtree = KDTree(X')
    idxs, dists = knn(kdtree, X', K+1, true)

    D = fill(Inf, (n, n))
    for i=1:n
        for (j, d) in zip(idxs[i], dists[i])
            D[i,j] = d
            D[j,i] = d
        end
    end
    D[diagind(D)] .= 0

    for k=1:n
        Dk = D[:,k] .+ D[k,:]'
        D[:,:] = min.(D, Dk)
    end
    return D
end

function ISOMAP(X, K)
    D = pathDistances(X, K)
    Y = MDS(D)
    return Y
end

function LLE(X, K)
    kdtree = KDTree(X')
    idxs, dists = knn(kdtree, X', K+1, true)

    tol = (K > 3) ? 1e-3 : 0

    W = zeros(n, n)
    for i=1:n
        js = idxs[i][2:end]
        Z = X[js,:] .- X[i,:]'
        C = Z * Z'
        C = C + I*(tr(C)*tol)
        ws = C \ ones(K)
        ws = ws / sum(ws)
        W[i,js] = ws
    end

    M = (W - I)' * (W - I)
    vals, vecs = eigs(M, nev=3, which=:SM)
    Y = vecs[:, 2:end]
    return Y
end


# Plot your embeddings
Y = ISOMAP(X, 12)
# Y = LLE(X, 12)
scene = scatter(Y[:,1], Y[:,2], color=X[:,3])

#### END HOMEWORK ASSIGNMENT ####
