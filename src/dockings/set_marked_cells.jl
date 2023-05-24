using BiochemicalAlgorithms
using Base.Threads

include("min_max_atoms.jl")

function set_marked_cells(atomballs::Vector{Meshes.Ball}, centroids::Array{Meshes.Point3,3}, roomcoordinates::Vector{Vector3{Float32}})
    # initialize vector with datatype of centroids
    colored_cells = Vector{Int64}()

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
    # store centroids that are inside a atom radius in colored_cells
    for i in CartesianIndices(centroids[min_x:max_x,min_y:max_y,min_z:max_z]) 
        for j in eachindex(atomballs)
            # move cartesian index i via min_x,min_y,min_z to get right index
            # I[] to get linear index of cartesian index
            index = I[CartesianIndex(min_x,min_y,min_z)+i]
            # check if centroid at index is in atomball j 
            if(Base.in(centroids[index],atomballs[j]))
                # stores indice of centroid if a centroid i lies
                # in an atomball j -> stored in colored_cells if not already in storage
                if(!Base.in(index, colored_cells))
                    push!(colored_cells,index)
                end
            end
        end
    end

    # return cells that represent protein in grid structure
    return colored_cells
end