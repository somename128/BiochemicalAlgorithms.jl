using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using ProgressBars
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("create_atomballs.jl")
include("correlation_docking.jl")
include("mass_center.jl")
include("rotate_atoms.jl")
include("helpers.jl")

N = Int32(64)
rotations = create_rotations()
R = Vector{Matrix3{Float32}}()
r = RotXYZ(deg2rad(80),deg2rad(0),deg2rad(0))
push!(R,r)
centroids = create_centroids(N, one(Int32))
protein_A = load_and_trans_pdb("dummy_protein_vol2.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
protein_B = load_and_trans_pdb("dummy_ligand_vol2.pdb", N)
roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
A = grid_representation(roomcoordiantes_atoms_A, N, centroids)
B = grid_representation(roomcoordiantes_atoms_B, N, centroids)

C = ifft(fft(A).*fft(B))

display(A)





