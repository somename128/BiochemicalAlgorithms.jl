using Distributed
using Base.Threads
using SharedArrays

include("min_max_atoms.jl")
include("transform_cells!.jl")

function set_grid(atomballs::Vector{Meshes.Ball{3,Float32}}, centroids::Array{Meshes.Point3f, 3}, roomcoordinates::Vector{Vector3{Float32}}, N::Int32, res::Int32, is_smaller::Bool)
    # initialize vector for storing index of colored cells
    grid = zeros(ComplexF32, N*res, N*res, N*res)
    # using shared array
    # grid = SharedArray{ComplexF32, 3}(Int64(N*res), Int64(N*res), Int64(N*res))

    #extract min max (in rounded int +/-2) of atom coordinates of protein
    min_max = min_max_atoms(roomcoordinates)
    min_x = min_max[1]
    max_x = min_max[2]
    min_y = min_max[3]
    max_y = min_max[4]
    min_z = min_max[5]
    max_z = min_max[6]

    # set grid cell where centroid inside atomball to one
    @threads for i in CartesianIndices(centroids[min_x*res:max_x*res,min_y*res:max_y*res,min_z*res:max_z*res])
        for j in eachindex(atomballs)
            # real index because indexing starts at one
            index = CartesianIndex(min_x*res,min_y*res,min_z*res)+i
            # check if centroid at index is in atomball j 
            if (Base.in(centroids[index],atomballs[j]))
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