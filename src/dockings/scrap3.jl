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
include("create_atomballs.jl")
include("load_trans_pdb.jl")
include("extract_roomcoordinates.jl")

# N = Int32(128)
# protein_A = load_and_trans_pdb("src/dockings/2ptc_protein.pdb", N)
# roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
inner_radius = Meshes.Ball(Meshes.Point(0,0,0),1)
outer_radius = Meshes.Ball(Meshes.Point(0,0,0),2)

p = Meshes.Point(0.9, 0.47, 0)

# !Base.in(p,inner_radius) && Base.in(p,outer_radius)
Base.in(p,inner_radius)