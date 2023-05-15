using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using ProgressBars
using Profile

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")



protein = load_and_trans_pdb("2ptc_ligand.pdb",128)

atomballs = create_atomballs(protein)
centroids = create_centroids(128,1)

@time colored_cells = set_marked_cells(atomballs,centroids,protein)