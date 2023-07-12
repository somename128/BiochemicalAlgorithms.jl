using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD2
using ProfileView
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations

include("quaternion_functions.jl")
include("create_rotations.jl")

rotations = create_rotations()

atoms_in_space_points = Vector{Vector3{Float32}}()
v = Vector3{Float32}(1,2,3)
push!(atoms_in_space_points, v)
w = Vector3{Float32}(3,4,5)
push!(atoms_in_space_points, w)

g(v::Vector3{Float32}) = rotate_vector(rotations[11], v)
atoms_rotated = g.(atoms_in_space_points)

typeof(atoms_rotated)

