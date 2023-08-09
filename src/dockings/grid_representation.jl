include("create_atomballs.jl")
include("set_marked_cells.jl")
include("create_inner_outer_grid.jl")
include("set_surface_cells.jl")
include("create_inner_outer_grid_vdW.jl")

function grid_representation(atoms::Vector{Tuple{String, Vector3{Float32}}}, gridsize::Int32, centroids::Array{Meshes.Point3f, 3}, res::Int32, is_smaller::Bool, vdW::Bool)

    # if vdW is enabled surface is defined as the difference "line" between two balls
    # called thinkness
    # inner points are defined as all points which are maximum [radius-thickness] 
    # angström away from atom point

    # extract roomcoordinates from tuple
    roomcoordinates = Vector{Vector3{Float32}}()
    [push!(roomcoordinates, i[2]) for i in atoms] 
    if (vdW)
        thickness = Float32(2)
        inner_radius = create_atomballs(atoms, -(thickness/2))
        outer_radius = create_atomballs(atoms, thickness/2)
        inner_cells = set_marked_cells(inner_radius, centroids, roomcoordinates, res)
        surface_cells = set_surface_cells(inner_radius, outer_radius, centroids, roomcoordinates, res)
        inner_outer_grid_3D = create_inner_outer_grid_vdW(inner_cells, surface_cells, gridsize, res)
    else
        # calculate atomballs around proteins 
        # zero because of no vdW surface calculation
        atomballs = create_atomballs(atoms, zero(Float32))
        # find all cells in grid which represent protein 
        colored_cells = set_marked_cells(atomballs, centroids, roomcoordinates, res)

        inner_outer_grid_3D = create_inner_outer_grid(colored_cells, gridsize, res, is_smaller)
    end

    return inner_outer_grid_3D
end