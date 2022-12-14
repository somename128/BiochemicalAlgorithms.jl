using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

grid = Meshes.CartesianGrid((0,0,0), (4,4,4), dims=(4,4,4))

viz(grid, showfacets = true)


