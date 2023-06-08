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

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("create_atomballs.jl")
include("correlation_docking.jl")

N = Int32(64)
protein = load_and_trans_pdb("dummy_ligand.pdb",N)
atoms = extract_roomcoordinates(protein)

centroids = create_centroids(N,one(Int32))


atomballs = create_atomballs(atoms)
colored_cells = set_marked_cells(atomballs, centroids, atoms)
grid = create_inner_outer_grid(colored_cells, N)

 # transfer atom coordinates in mesh points
 atoms_in_space_points = Base.Vector{Meshes.Point3}()

 for i in CartesianIndices(grid)
    if(grid[i] != 0)
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end
end

min_max_atoms(atoms)



