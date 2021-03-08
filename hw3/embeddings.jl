using DataFrames
using CSV
using Makie

include("utils.jl")

filename = "data/swissroll.txt"
data = CSV.read(filename, DataFrame; datarow=1, delim=' ', ignorerepeated=true)
X = [data.Column1 data.Column2 data.Column3]

scene = scatter(X[:,1], X[:,2], X[:,3], color=X[:,3], markersize=0.3)

#### ADD CODE FOR YOUR TWO METHODS HERE ####

# Plot your embeddings

#### END HOMEWORK ASSIGNMENT ####
