include("utils.jl")

filename = "meshes/moomoo.off"
# filename = "meshes/166.off"
X, T = readoff(filename)
nv = size(X, 1)
nt = size(T, 1)

# FILL IN FUNCTIONS TO COMPUTE GAUSSIAN CURVATURE HERE ##########
# Feel free to rename these functions
function gaussianCurvature1(X, T)
    return rand(nv)
end

function gaussianCurvature2(X, T)
    return rand(nv)
end
# END HOMEWORK PROBLEM ######

scene1 = showdescriptor(X, T, gaussianCurvature1(X, T))
scene2 = showdescriptor(X, T, gaussianCurvature2(X, T))
