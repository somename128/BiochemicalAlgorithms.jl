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



protein = load_and_trans_pdb("2ptc_protein.pdb")

atomballs = create_atomballs(protein)
centroids = create_centroids(128,1)
#=
Profile.clear()
@profile set_marked_cells(atomballs,centroids,protein)
ProfileView.view()
=#


min_max = min_max_atoms(protein)

min_x = min_max[1]
max_x = min_max[2]
min_y = min_max[3]
max_y = min_max[4]
min_z = min_max[5]
max_z = min_max[6]

colored_cells_1 = Vector{Int64}()
colored_cells_2 = Vector{Int64}()
#=
@time begin
colored_cells_1 = set_marked_cells(atomballs,centroids,protein)
end
=#
I = LinearIndices(centroids)

@time begin
# store centroids that are inside a atom radius in colored_cells
colored_cells_1 = set_marked_cells(atomballs,centroids,protein)
end

length(colored_cells_1)

#=
println(length(colored_cells_1))
println(length(colored_cells_2))
C = CartesianIndices(centroids)

println(colored_cells_2[100])
println(I[colored_cells_2[100]])
println(C[colored_cells_2[100]])

for i in CartesianIndices(centroids[min_x:max_x,min_y:max_y,min_z:max_z])
    index = I[(i+CartesianIndex(min_x,min_y,min_z))]
    println(i+CartesianIndex(min_x,min_y,min_z))
    println(I[(i+CartesianIndex(min_x,min_y,min_z))])
    println(centroids[index])
end
=#

