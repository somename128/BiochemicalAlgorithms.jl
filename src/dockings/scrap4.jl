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



protein = load_and_trans_pdb("2ptc_ligand.pdb")

atomballs = create_atomballs(protein)
centroids = create_centroids(128,1)
#=
Profile.clear()
@profile set_marked_cells(atomballs,centroids,protein)
ProfileView.view()
=#

#extract min max (in rounded int) of atom coordinates of protein
min_max = min_max_atoms(protein)
min_x = min_max[1]
max_x = min_max[2]
min_y = min_max[3]
max_y = min_max[4]
min_z = min_max[5]
max_z = min_max[6]

colored_cells = Vector{Int64}()

# store centroids that are inside a atom radius in colored_cells
@time begin
for i in eachindex(centroids), j in eachindex(atomballs)
    if(Base.in(centroids[i],atomballs[j]))
        # returns vector thats why position[1]
        # dont know if vector of vectors or number better 
        # for future calculations
        #
        # findall returns indice of i in centroids if a centroid i lies
        # in an atomball j -> stored in colored_cells if not already in storage
        position = i
        println(position)
        if(!Base.in(position, colored_cells))
            push!(colored_cells,position[1])
        end
    end
end
end

length(colored_cells)
