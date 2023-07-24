include("min_max_atoms.jl")

function set_surface_cells(inner_radius::Vector{Meshes.Ball}, outer_radius::Vector{Meshes.Ball}, centroids::Array{Meshes.Point3f,3}, roomcoordinates::Vector{Vector3{Float32}}, res::Int32)
    # initalize vector for storing index of surface cells
    surface_cells = Vector{Int32}()

    #extract min max (in rounded int +/-2) of atom coordinates of protein
    min_max = min_max_atoms(roomcoordinates)
    min_x = min_max[1]
    max_x = min_max[2]
    min_y = min_max[3]
    max_y = min_max[4]
    min_z = min_max[5]
    max_z = min_max[6]

    # extract LinearIndices from centroids
    # still not sure how this stuff works
    I = LinearIndices(centroids)
    # println("Set marked cells...")
    # store centroids that are not inside the inner radius but inside the outer radius
    # in between of the two radians therefore on the surface
    # inner and outer radius should be of the same length because of the constrution process
    for i in CartesianIndices(centroids[min_x*res:max_x*res,min_y*res:max_y*res,min_z*res:max_z*res]), j in eachindex(inner_radius)
        # move cartesian index i via min_x,min_y,min_z to get right index
        # I[] to get linear index of cartesian index
        index = I[CartesianIndex(min_x*res,min_y*res,min_z*res)+i]
        # check if centroid at index is in atomball j 
        if(!Base.in(centroids[index],inner_radius[j]) && Base.in(centroids[index],outer_radius[j]))
            # stores indice of centroid if a centroid i lies
            # in the surface area j -> stored in surface_cells if not already in storage
            if(!Base.in(index, surface_cells))
                push!(surface_cells, index)
            end
        end
    end

    # return cells that represent protein in grid structure
    return surface_cells
end