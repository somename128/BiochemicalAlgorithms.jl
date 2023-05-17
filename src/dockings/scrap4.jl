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
include("generate_record.jl")

N = 128
protein = load_and_trans_pdb("2ptc_protein.pdb",N)
centroids = create_centroids(N,1)
A = grid_representation(protein,N, centroids)
rotations = create_rotations()


# Profile.clear()
@time record = generate_record(A, rotations[100], "2ptc_ligand.pdb", centroids,N)
# ProfileView.view()
N

