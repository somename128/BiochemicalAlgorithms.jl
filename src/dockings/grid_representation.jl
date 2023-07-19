include("create_atomballs.jl")
include("set_marked_cells.jl")
include("create_inner_outer_grid.jl")
include("set_surface_cells.jl")
include("create_inner_outer_grid_vdW.jl")

function grid_representation(atoms::Vector{Vector3{Float32}}, gridsize::Int32, centroids::Array{Meshes.Point3, 3}, is_smaller::Bool, vdW::Bool)
    
    if (vdW)
        # print("You choose the vdW-surface option. How thick should the surface be? ")
        # thickness = readline()
        thickness = Float32(1)
        inner_radius = create_atomballs(atoms, Float32(1.8))
        outer_radius = create_atomballs(atoms, Float32(1.8) + thickness)
        inner_cells = set_marked_cells(inner_radius, centroids, atoms)
        surface_cells = set_surface_cells(inner_radius, outer_radius, centroids, atoms)
        inner_outer_grid_3D = create_inner_outer_grid_vdW(inner_cells, surface_cells, gridsize)
    else
        # calculate atomballs around proteins atoms with
        # radius r
        atomballs = create_atomballs(atoms, Float32(0.7))
        # find all cells in grid which represent protein 
        colored_cells = set_marked_cells(atomballs, centroids, atoms)

        inner_outer_grid_3D = create_inner_outer_grid(colored_cells, gridsize, is_smaller)
    end

    return inner_outer_grid_3D
end