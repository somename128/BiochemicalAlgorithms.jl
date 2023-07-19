using Meshes

function create_centroids(gridsize::Int32, resolution::Int32)
    # initialize vector for points 
    # size is gridsize*resolution and
    # we would like to have only the centroids of the cells
    # not he "boundary" points
    centroids = Vector{Meshes.Point3f}()
    resfactor = 1/resolution

    for k in 1:gridsize*resolution, j in 1:gridsize*resolution, i in 1:gridsize*resolution
        c = Meshes.Point3f((i-1)*resfactor + resfactor/2, (j-1)*resfactor + resfactor/2, 
            (k-1)*resfactor + resfactor/2)
        push!(centroids, c)
    end

    return centroids
end