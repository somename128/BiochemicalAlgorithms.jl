using Meshes, MeshViz
using Makie, WGLMakie
using ProfileView
using Profile
using BenchmarkTools
using JET
using BiochemicalAlgorithms

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("create_rotations.jl")
include("create_centroids.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("set_gridsize.jl")
include("create_rotations2.jl")
include("create_atomballs.jl")
include("set_surface_cells_fast.jl")
include("set_surface_cells.jl")
include("set_marked_cells_fast.jl")
include("set_marked_cells.jl")
include("rotate_atoms.jl")
include("create_inner_outer_grid.jl")


pathA = "src/dockings/testproteins/2ptc_proteinase.pdb"
pathB = "src/dockings/testproteins/2ptc_inhibitor.pdb"
# init
res = Int32(1)
N = set_gridsize(pathA, pathB)
println("Gridsize: ", N)
# N = Int32(32)
# load and translate protein a
protein_A = load_and_trans_pdb(pathA, N)
atoms = extract_roomcoordinates(protein_A)
# calculate centroids in a NxNxN grid with cells
# of 1 angström, only done once
centroids = create_centroids(N, res)
thickness = Float32(2)
inner_radius = create_atomballs(atoms, -(thickness/2))
outer_radius = create_atomballs(atoms, thickness/2)
roomcoordinates = Array{Vector3{Float32}}(undef, length(atoms))
[roomcoordinates[i] = atoms[i][2] for i in eachindex(atoms)]
rotations = create_rotations()


# @profview set_surface_cells_fast(inner_radius, outer_radius, centroids, roomcoordinates, res)

# @time inner_fast = set_marked_cells_fast(inner_radius, centroids, roomcoordinates, res)
# @time surf = set_surface_cells(inner_radius, outer_radius, centroids, roomcoordinates, res)
inner = set_marked_cells(inner_radius, centroids, roomcoordinates, res)

# println(sort(surf) == sort(surf_fast))
# println(sort(inner) == sort(inner_fast))
@btime create_inner_outer_grid(inner, N, res, false)
