include("create_atomballs.jl")
include("set_marked_cells.jl")
include("create_inner_outer_grid.jl")

function grid_representation(atoms::Vector{Vector3{Float32}}, gridsize::Int32, centroids::Array{Meshes.Point3, 3}, is_smaller::Bool)
    # calculate atomballs around proteins atoms
    atomballs = create_atomballs(atoms)
    # find all cells in grid which represent protein 
    colored_cells = set_marked_cells(atomballs,centroids,atoms)

    inner_outer_grid_3D = create_inner_outer_grid(colored_cells, gridsize, is_smaller)

    return inner_outer_grid_3D
end