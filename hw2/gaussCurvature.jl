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

function gaussianCurvature2(X, T)
    return rand(nv)
end
# END HOMEWORK PROBLEM ######

scene1 = showdescriptor(X, T, gaussianCurvature1(X, T))
scene2 = showdescriptor(X, T, gaussianCurvature2(X, T))
