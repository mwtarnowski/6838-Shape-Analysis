module ElasticRods

export elasticRods

using LinearAlgebra
using SparseArrays
using Rotations
using WGLMakie

function multicross(x, y)
    reduce(hcat, cross.(eachcol(x), eachcol(y)))
end

mutable struct CurveData
    verts :: Array{Float64}
    velocities :: Array{Float64}
    totalTwist :: Float64
    edgesR :: Array{Float64}
    edgesL :: Array{Float64}
    midpointsR :: Array{Float64}
    edgeLengthsR :: Array{Float64}
    edgeLengthsL :: Array{Float64}
    dualLengths :: Array{Float64}
    massMatrix :: Array{Float64}
    totalLength :: Float64
    tangentsR :: Array{Float64}
    tangentsL :: Array{Float64}
    curvatureBinormalsDenom :: Array{Float64}
    curvatureBinormals :: Array{Float64}
    curvatureSquared :: Array{Float64}
    parallelTransport :: Array{RotMatrix{3, Float64}}
    u0 :: Array{Float64}
    v0 :: Array{Float64}
    u :: Array{Float64}
    v :: Array{Float64}
    bendingDensity :: Array{Float64}
    bendingEnergy :: Float64
    twistEnergy :: Float64

    function CurveData(verts, velocities, totalTwist)
        curveData = new(verts, velocities, totalTwist)

        curveData.edgesR = circshift(curveData.verts, (0, -1)) .- curveData.verts
        curveData.edgesL = circshift(curveData.edgesR, (0, 1))
        curveData.midpointsR = 0.5 .* (circshift(curveData.verts, (0, -1)) .+ curveData.verts)
        curveData.edgeLengthsR = sqrt.(sum(curveData.edgesR.^2, dims=1))
        curveData.edgeLengthsL = circshift(curveData.edgeLengthsR, (0, 1))
        curveData.dualLengths = 0.5 .* (curveData.edgeLengthsR .+ curveData.edgeLengthsL)
        curveData.totalLength = sum(curveData.dualLengths)
        curveData.massMatrix = kron(Diagonal(curveData.dualLengths[:]), Diagonal(ones(3)))
        
        curveData.tangentsR = curveData.edgesR ./ curveData.edgeLengthsR
        curveData.tangentsL = circshift(curveData.tangentsR, (0, 1))
        curveData.curvatureBinormalsDenom = curveData.edgeLengthsL .* curveData.edgeLengthsR .+
                                            sum(curveData.edgesL .* curveData.edgesR, dims=1)
        curveData.curvatureBinormals = 2 * multicross(curveData.edgesL, curveData.edgesR) ./
                                       curveData.curvatureBinormalsDenom
        curveData.curvatureSquared = sum(curveData.curvatureBinormals.^2, dims=1)
        
        ### PROBLEM 3(a) - YOUR CODE HERE
        curveData.parallelTransport = fill(one(RotMatrix{3, Float64}), size(verts, 2))
        ### END HOMEWORK PROBLEM
        
        curveData.u0 = zeros(size(verts))
        curveData.v0 = zeros(size(verts))
        curveData.u0[:, 1] = cross(curveData.tangentsR[:, 1], [0;0;1])
        curveData.v0[:, 1] = cross(curveData.u0[:, 1], curveData.tangentsR[:, 1])
        for j = 2:size(curveData.u0, 2)
            P = curveData.parallelTransport[j]
            curveData.u0[:, j] = P * curveData.u0[:, j - 1]
            curveData.v0[:, j] = P * curveData.v0[:, j - 1]
        end
        cumTwist = curveData.totalTwist * cumsum([0 curveData.dualLengths[2:end]'], dims=2) ./ curveData.totalLength
        curveData.u = cos.(cumTwist) .* curveData.u0 .+ sin.(cumTwist) .* curveData.v0
        curveData.v = cos.(cumTwist) .* curveData.v0 .- sin.(cumTwist) .* curveData.u0
                                    
        curveData.bendingDensity = sum(curveData.curvatureBinormals.^2, dims=1) ./ curveData.dualLengths
        curveData.bendingEnergy = sum(curveData.bendingDensity)
        curveData.twistEnergy = curveData.totalTwist.^2 ./ curveData.totalLength
        return curveData
    end
end

function elasticRods(bendModulus = 1, twistModulus = 1, totalTwist = pi)
    nSamples = 100
    dt = 0.001
    curveFunction(t) = [cos.(2pi .* t); sin.(2pi .* t); 0.3 .* sin.(4pi .* t)]
    verts = curveFunction(range(0, stop=1, length=nSamples + 1)')[:, 1:nSamples]
    curveData = CurveData(verts, zeros(3, nSamples), totalTwist)

    links = [(1:nSamples) circshift((1:nSamples), -1)]
    DCii = repeat((1:nSamples), 1, 6)
    DCjj = mod.(3 .* ((1:nSamples) .- 1) .+ (0:5)', 3 .* nSamples) .+ 1

    # Observables to be plotted
    loop = Node(Point3.(eachcol([verts verts[:, 1]])))
    bases = Node(Point3.(eachcol(curveData.midpointsR)))
    vecs = Node(Point3.(eachcol(curveData.u)))

    fig = Figure(resolution=(1000, 1000))
    ax = LScene(fig)
    lines!(loop, linewidth=30)
    arrows!(bases, vecs, linewidth=10, linecolor=:red, arrowcolor=:red, arrowsize=0.01, lengthscale=0.1)
    oax = ax.scene[OldAxis]
    oax.showgrid = (false, false, false)
    oax.showticks = (false, false, false)
    fig[1, 1] = ax
    display(fig)

    for i = 1:10000
        bendForce = computeBendForce(curveData)
        twistForce = computeTwistForce(curveData)
        totalForce = bendModulus .* bendForce .+ twistModulus .* twistForce

        verts, velocities = symplecticEuler(curveData.verts, curveData.velocities, totalForce, dt)
        verts, velocities = fastProjection(curveData.verts, verts, curveData.massMatrix, curveData.edgeLengthsR, dt, DCii, DCjj)
        curveData = CurveData(verts, velocities, totalTwist)
        
        if mod(i - 1, 10) == 0
            sleep(0.001)
            loop[] = Point3.(eachcol([verts verts[:, 1]]))
            bases[] = Point3.(eachcol(curveData.midpointsR))
            vecs[] = Point3.(eachcol(curveData.u))
        end
    end
end

### PROBLEM 3(c) Part I - YOUR CODE HERE
function computeBendForce(curveData)
    bendForce = zeros(size(curveData.verts));
    return bendForce
end
### END HOMEWORK PROBLEM

### PROBLEM 3(c) Part II - YOUR CODE HERE
function computeTwistForce(curveData)
    twistForce = zeros(size(curveData.verts));
    return twistForce
end
### END HOMEWORK PROBLEM

function symplecticEuler(verts0, velocities0, totalForce, dt)
    # Update velocities and positions via Symplectic Euler Integration
    velocities = velocities0 + dt * totalForce;
    verts = verts0 + dt * velocities;
    return verts, velocities
end

function fastProjection(verts0, verts, massMatrix, lengthsR, dt, DCii, DCjj)
    # Project positions to satisfy length constraints
    edgesR = circshift(verts, (0, -1)) - verts;
    constraint = (sum(edgesR.^2, dims=1) .- lengthsR.^2)';
    while maximum(abs.(constraint)) > 1e-10
        constraintGrad = 2 * sparse(DCii[:], DCjj[:], vec([-edgesR' edgesR']), length(constraint), length(verts));
        MinvDC = massMatrix \ constraintGrad';
        DCMinvDC = constraintGrad * MinvDC;
        dLambda = DCMinvDC \ constraint;
        dx = -reshape(MinvDC * dLambda, size(verts));
        verts = verts + dx;
        edgesR = circshift(verts, (0, -1)) .- verts;
        constraint = (sum(edgesR.^2, dims=1) .- lengthsR.^2)';
    end
    velocities = (verts .- verts0) ./ dt;
    return verts, velocities
end

end
