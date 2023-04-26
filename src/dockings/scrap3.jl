using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

test = Vector{Int32}([1,2,3,4])

insert = Vector{Int32}([1,2,5,6,7,4,3,2,1,2])

for (index, value) in insert
    println(index,"",value)
end




