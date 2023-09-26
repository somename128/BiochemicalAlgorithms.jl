using Distributed
using Base.Threads
using SharedArrays
using Meshes

include("transform_cells!.jl")

function set_grid(atomballs::Vector{Meshes.Ball{3,Float32}}, centroids::Array{Meshes.Point3f, 3}, roomcoordinates::Vector{Vector3{Float32}}, N::Int32, res::Int32, is_smaller::Bool)
    # initialize vector for storing index of colored cells
    grid = zeros(ComplexF32, N*res, N*res, N*res)
    # using shared array
    # grid = SharedArray{ComplexF32, 3}(Int64(N*res), Int64(N*res), Int64(N*res))

    # set grid cell where centroid inside atomball to one
    @threads for i in eachindex(atomballs)
        min_x = Int32(round(atomballs[i].center.coords[1])) - Int32(3)
        max_x = Int32(round(atomballs[i].center.coords[1])) + Int32(3)
        min_y = Int32(round(atomballs[i].center.coords[2])) - Int32(3)
        max_y = Int32(round(atomballs[i].center.coords[2])) + Int32(3)
        min_z = Int32(round(atomballs[i].center.coords[3])) - Int32(3)
        max_z = Int32(round(atomballs[i].center.coords[3])) + Int32(3)

        for j in CartesianIndices(centroids[min_x*res:max_x*res,min_y*res:max_y*res,min_z*res:max_z*res])
            # real index because indexing starts at one
            index = CartesianIndex(min_x*res,min_y*res,min_z*res)+j
            # check if centroid at index is in atomball j 
            if (Meshes.in(centroids[index],atomballs[i]))
                grid[index] = one(ComplexF32)
            end
        end
    end
    
    # return grid structure
    if (!is_smaller)
        return transform_cells!(grid)
    else
        return grid
    end
end