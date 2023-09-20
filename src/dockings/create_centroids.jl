using Meshes

function create_centroids(gridsize::Int32, res::Int32)
    # initialize array for points 
    # size is gridsize*resolution and
    # we would like to have only the centroids of the cells
    # not the "boundary" points
    M = gridsize * res
    centroids = Array{Meshes.Point3f}(undef, M, M, M)
    resfactor = 1/res

    for i in CartesianIndices(centroids)
        centroids[i] = Meshes.Point3f((i[1]-1)*resfactor + resfactor/2, (i[2]-1)*resfactor + resfactor/2, 
            (i[3]-1)*resfactor + resfactor/2)
    end

    return centroids
end