using Arpack
using Printf
using SparseArrays

include("utils.jl")

filename = "data/moomoo.off"
# filename = "data/166.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)


triangleArea(X1, X2, X3) = norm(cross(X1-X2, X3-X2)) / 2

# ADD CODE TO COMPUTE SURFACE AREA HERE #############
function surfaceArea(X, T)
    area = 0
    for i=1:nt
        v1,v2,v3 = T[i,:]
        X1 = X[v1,:]
        X2 = X[v2,:]
        X3 = X[v3,:]
        A = triangleArea(X1, X2, X3)
        area += A
    end
    return area
end
# END HOMEWORK ASSIGNMENT #########

@printf("The surface area of %s is %f\n", filename, surfaceArea(X, T))

# ADD CODE TO COMPUTE COTANGENT LAPLACIAN HERE ##########
function cotLaplacian(X, T)
    cots = zeros(nt, 3)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        e12 = X[v2,:] - X[v1,:]
        e23 = X[v3,:] - X[v2,:]
        e31 = X[v1,:] - X[v3,:]
        cots[i,1] = -dot(e12, e31) / norm(cross(e12, e31))
        cots[i,2] = -dot(e23, e12) / norm(cross(e23, e12))
        cots[i,3] = -dot(e31, e23) / norm(cross(e31, e23))
    end

    I = [T[:,1]; T[:,2]; T[:,3]]
    J = [T[:,2]; T[:,3]; T[:,1]]
    V = [-cots[:,3]; -cots[:,1]; -cots[:,2]]
    L = sparse(I, J, V, nv, nv)
    L = L + L'

    Vdiag = vec(sum(L, dims=2))
    L = L - spdiagm(0 => Vdiag)
    return L
end
# END HOMEWORK ASSIGNMENT #######

# Sanity checks: Laplacian is symmetric and positive definite
L = cotLaplacian(X, T)
eigenvals, _ = eigs(L, nev=10, which=:SM)
println(eigenvals)
println(norm(L - L'))

# ADD CODE FOR DIVIDED DIFFERENCES HERE ######
function dividedDifferences(X, T)
    gradApprox = zeros(nv, 3)
    h = 0.001
    hI = h * Matrix(I, 3, 3)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        X1 = X[v1,:]
        X2 = X[v2,:]
        X3 = X[v3,:]
        for k=1:3
            hv = hI[k,:]
            gradApprox[v1,k] += (triangleArea(X1+hv, X2,    X3   ) - triangleArea(X1-hv, X2,    X3   )) / (2*h)
            gradApprox[v2,k] += (triangleArea(X1,    X2+hv, X3   ) - triangleArea(X1,    X2-hv, X3   )) / (2*h)
            gradApprox[v3,k] += (triangleArea(X1,    X2,    X3+hv) - triangleArea(X1,    X2,    X3-hv)) / (2*h)
        end
    end
    return gradApprox
end
# END HOMEWORK ASSIGNMENT #########

# Check that gradApprox and .5*L*X are similar
println(norm(.5 .* (L*X) - dividedDifferences(X,T)))

# ADD CODE FOR COMPUTING THE BARYCENTRIC AREA VECTOR HERE ######
function barycentricArea(X, T)
    M = zeros(nv)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        X1 = X[v1,:]
        X2 = X[v2,:]
        X3 = X[v3,:]
        A = triangleArea(X1, X2, X3)
        M[v1] += A / 3
        M[v2] += A / 3
        M[v3] += A / 3
    end
    return M
end
# END HOMEWORK ASSIGNMENT ########

# ADD CODE FOR COMPUTING POINTWISE MEAN CURVATURE HERE #####
function meanCurvature(X, T)
    L = cotLaplacian(X, T)
    M = barycentricArea(X, T)
    Hn = 0.5 * (L * X) ./ M

    normal = zeros(nv, 3)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        e12 = X[v2,:] - X[v1,:]
        e23 = X[v3,:] - X[v2,:]
        e31 = X[v1,:] - X[v3,:]
        n = cross(e12, e23)  # weighted by the triangle area
        normal[v1,:] += n
        normal[v2,:] += n
        normal[v3,:] += n
    end
    normal = mapslices(normalize, normal, dims=[2])

    H = zeros(nv)
    for i=1:nv
        H[i] = dot(Hn[i,:], normal[i,:])
    end
    return H
end
# END HOMEWORK ASSIGNMENT #######

H = meanCurvature(X, T)
scene = showdescriptor(X, T, H)


## Mean curvature flow ##
function curvatureFlowEuler(X, T)
    Xt = copy(X)
    maxiters = 1000
    # ADD CODE FOR THE EXPLICIT INTEGRATOR HERE ####
    tau = 5e-7 * surfaceArea(Xt, T)
    L = cotLaplacian(Xt, T)
    for t=1:maxiters
        # L = cotLaplacian(Xt, T)
        M = barycentricArea(Xt, T)
        Xt[:] = Xt - tau * (L * Xt) ./ M
    end
    # END HOMEWORK ASSIGNMENT #####
    H = meanCurvature(Xt, T)
    # Uncomment to show mean curvature at the end
    # scene = showdescriptor(Xt, T, H)
end

function curvatureFlowImplicit(X, T)
    Xt = copy(X)
    maxiters = 1000
    # ADD CODE FOR SEMI-IMPLICIT INTEGRATOR HERE ####
    tau = 5e-7 * surfaceArea(Xt, T)
    L = cotLaplacian(Xt, T)
    for t=1:maxiters
        # L = cotLaplacian(Xt, T)
        M = barycentricArea(Xt, T)
        Xt[:] = (spdiagm(0 => M) + tau * L) \ (M .* Xt)
    end
    # END HOMEWORK ASSIGNMENT
    H = meanCurvature(Xt, T)
    # Uncomment to show mean curvature at the end
    # scene = showdescriptor(Xt, T, H)
end

# Uncomment either the implicit or explicit flow to see your results
# curvatureFlowEuler(X,T)
# curvatureFlowImplicit(X,T)
