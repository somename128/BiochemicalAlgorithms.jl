using BiochemicalAlgorithms

include("create_atomballs.jl")
include("create_centroids.jl")
include("set_marked_cells.jl")
include("create_inner_outer_grid.jl")

function grid_representation(protein::Molecule{Float32}, gridsize::Int64, centroids::Array{Meshes.Point3, 3})
    # N for grid size
    N = gridsize
    # calculate atomballs around proteins atoms
    atomballs = create_atomballs(protein)
    # find all cells in grid which represent protein 
    colored_cells = set_marked_cells(atomballs,centroids,protein)

    inner_outer_grid_3D = create_inner_outer_grid(colored_cells, gridsize)

    return inner_outer_grid_3D
end