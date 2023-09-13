using Meshes, MeshViz
using Makie, WGLMakie
using ProfileView
using Profile
using BenchmarkTools

include("extract_max.jl")
include("extract_max_fast.jl")
include("create_rotations.jl")
include("create_rotations2.jl")

gridsize = 1
res = 2
B = Array{ComplexF32,3}(undef, gridsize*res, gridsize*res, gridsize*res)