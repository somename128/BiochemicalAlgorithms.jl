using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

test = Vector{Int32}([1,2,3,4])

insert = Vector{Int32}([1,2,5,6,7,4,3,2,1,2])

for i in insert
    if(!Base.in(i,test))
        push!(test,i)
    end
end

test

