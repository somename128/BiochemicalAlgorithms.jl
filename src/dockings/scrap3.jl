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
include("get_degrees.jl")

rotations = create_rotations()

atoms_in_space_points = Vector{Vector3{Float32}}()
v = Vector3{Float32}(1,2,3)
push!(atoms_in_space_points, v)
w = Vector3{Float32}(3,4,5)
push!(atoms_in_space_points, w)

# g(v::Vector3{Float32}) = rotate_vector(rotations[11], v)
# atoms_rotated = g.(atoms_in_space_points)
# α degrees rotation aroud x-axis
q1 = quat_from_axisangle([1,0,0], deg2rad(40))
# β degrees rotation aroud y-axis
q2 = quat_from_axisangle([0,1,0], deg2rad(20))
# β degrees rotation aroud z-axis
q3 = quat_from_axisangle([0,0,1], deg2rad(20))
# combine quaternions for full rotation around xyz
q = q1*q2*q3

R = rotmatrix_from_quat(q)
typeof(R)

