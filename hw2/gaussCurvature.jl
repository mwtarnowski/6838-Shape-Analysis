include("utils.jl")

filename = "data/moomoo.off"
# filename = "data/166.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)

# FILL IN FUNCTIONS TO COMPUTE GAUSSIAN CURVATURE HERE ##########
function voronoiRegionArea(X, T)
    area = zeros(nv)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        e12 = X[v2,:] - X[v1,:]
        e23 = X[v3,:] - X[v2,:]
        e31 = X[v1,:] - X[v3,:]
        cot1 = -dot(e12, e31) / norm(cross(e12, e31))
        cot2 = -dot(e23, e12) / norm(cross(e23, e12))
        cot3 = -dot(e31, e23) / norm(cross(e31, e23))
        area[v1] += (dot(e31, e31)*cot2 + dot(e12, e12)*cot3) / 8
        area[v2] += (dot(e12, e12)*cot3 + dot(e23, e23)*cot1) / 8
        area[v3] += (dot(e23, e23)*cot1 + dot(e31, e31)*cot2) / 8
    end
    return area
end

function barycentricRegionArea(X, T)
    area = zeros(nv)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        e12 = X[v2,:] - X[v1,:]
        e23 = X[v3,:] - X[v2,:]
        e31 = X[v1,:] - X[v3,:]
        A = norm(cross(e12, e23)) / 2
        area[v1] += A / 3
        area[v2] += A / 3
        area[v3] += A / 3
    end
    return area
end

function vertexNormal(X, T)
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
    return normal
end

# M.Meyer, M.Desbrun, P.Schroder, A.H.Barr,
# Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.
function gaussianCurvature1(X, T)
    defect = fill(2*pi, nv)
    for i=1:nt
        v1,v2,v3 = T[i,:]
        e12 = X[v2,:] - X[v1,:]
        e23 = X[v3,:] - X[v2,:]
        e31 = X[v1,:] - X[v3,:]
        defect[v1] -= acos(-dot(e12, e31) / (norm(e12) * norm(e31)))
        defect[v2] -= acos(-dot(e23, e12) / (norm(e23) * norm(e12)))
        defect[v3] -= acos(-dot(e31, e23) / (norm(e31) * norm(e23)))
    end
    # curvature = defect ./ barycentricRegionArea(X, T)
    curvature = defect ./ voronoiRegionArea(X, T)
    return curvature
end

# G.Taubin,
# Estimating the Tensor of Curvature of a Surface from a Polyhedral Approximation.
function gaussianCurvature2(X, T)
    normal = vertexNormal(X, T)

    M = zeros(nv, 3, 3)
    ws = zeros(nv)
    for i = 1:nt
        v1,v2,v3 = T[i,:]
        e12 = X[v2,:] - X[v1,:]
        e23 = X[v3,:] - X[v2,:]
        e31 = X[v1,:] - X[v3,:]
        A = norm(cross(e12, e23)) / 2
        n1 = normal[v1,:]
        n2 = normal[v2,:]
        n3 = normal[v3,:]
        # vertex 1
        P1 = I - n1*n1'
        t12 = -normalize(P1 * e12)
        t13 =  normalize(P1 * e31)
        k12 =  2*dot(n1, e12) / dot(e12, e12)
        k13 = -2*dot(n1, e31) / dot(e31, e31)
        M[v1,:,:] += (k12 * (t12*t12') + k13 * (t13*t13')) * A
        ws[v1] += 2*A
        # vertex 2
        P2 = I - n2*n2'
        t23 = -normalize(P2 * e23)
        t21 =  normalize(P2 * e12)
        k23 =  2*dot(n2, e23) / dot(e23, e23)
        k21 = -2*dot(n2, e12) / dot(e12, e12)
        M[v2,:,:] += (k21 * (t21*t21') + k23 * (t23*t23')) * A
        ws[v2] += 2*A
        # vertex 3
        P3 = I - n3*n3'
        t31 = -normalize(P3 * e31)
        t32 =  normalize(P3 * e23)
        k31 =  2*dot(n3, e31) / dot(e31, e31)
        k32 = -2*dot(n3, e23) / dot(e23, e23)
        M[v3,:,:] += (k31 * (t31*t31') + k32 * (t32*t32')) * A
        ws[v3] += 2*A
    end
    M = M ./ ws

    curvature = zeros(nv)
    for i=1:nv
        l = eigvals(M[i,:,:], sortby = x -> -abs(x))
        curvature[i] = (3*l[1] - l[2]) * (3*l[2] - l[1])
    end
    return curvature
end
# END HOMEWORK PROBLEM ######

scene1 = showdescriptor(X, T, gaussianCurvature1(X, T))
scene2 = showdescriptor(X, T, gaussianCurvature2(X, T))
