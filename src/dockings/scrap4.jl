using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using JLD

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")


protein = load_and_trans_pdb("2ptc_protein.pdb")

atomballs = create_atomballs(protein)
centroids = create_centroids(128,1)

#extract min max (in rounded int) of atom coordinates of protein
min_max = min_max_atoms(protein)
min_x = min_max[1]
max_x = min_max[2]
min_y = min_max[3]
max_y = min_max[4]
min_z = min_max[5]
max_z = min_max[6]

colored_cells = Vector{Meshes.CartesianIndex{3}}()

for i in centroids[min_x:max_x,min_y:max_y,min_z:max_z], j in atomballs
    if(Base.in(i,j))
        # returns vector thats why position[1]
        # dont know if vector of vectors or number better 
        # for future calculations
        #
        # findall returns indice of i in centroids if a centroid i lies
        # in an atomball j -> stored in colored_cells if not already in storage
        position = findall(item -> item == i, centroids)
        if(!Base.in(position[1], colored_cells))
            push!(colored_cells,position[1])
        end
    end  
end
colored_cells
#=
# grid properties
lower_left = (0,0,0)
N = 128
spcing_factor = 1
upper_right = (N,N,N) 
spcing = (N*spcing_factor, N*spcing_factor, N*spcing_factor)

# create 3D grid with N/N*spcing_factor spacing in each dimension, origin at (0,0,0)
grid = Meshes.CartesianGrid(lower_left, upper_right, dims=spcing)
# centroid of each cell
println("Build centroids...")
centroids = Meshes.centroid.(grid)

centroids = reshape([centroids...], N, N, N)

for i in centroids[min_x:max_x,min_y:max_y,min_z:max_z]
    position = findall(item -> item == i, centroids)
    println(position[1])
end

# viz(centroids[1])
=#

