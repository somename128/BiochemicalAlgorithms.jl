using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using ProgressBars
using Profile
using TypedTables

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("create_atomballs.jl")
include("correlation_docking.jl")

N = 128
protein = load_and_trans_pdb("2ptc_ligand.pdb",N)
atoms = extract_roomcoordinates(protein)
centroids = create_centroids(N,1)

@time begin
atomballs = create_atomballs(atoms)
colored_cells = set_marked_cells(atomballs, centroids, atoms)
grid = create_inner_outer_grid(colored_cells, N)
end

println("yey")

