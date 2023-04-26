using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")



# protein = load_and_trans_pdb("2ptc_ligand.pdb")

# atomballs = create_atomballs(protein)
# println(length(atomballs))
centroids = create_centroids(128,1)
println(length(centroids))

min_max = min_max_atoms(protein)
    min_x = min_max[1]
    max_x = min_max[2]
    min_y = min_max[3]
    max_y = min_max[4]
    min_z = min_max[5]
    max_z = min_max[6]

relevant_centroids = Vector{Meshes.Point3}()

for i in centroids[min_x:max_x,min_y:max_y,min_z:max_z]
    push!(relevant_centroids,i)
end
length(relevant_centroids)


