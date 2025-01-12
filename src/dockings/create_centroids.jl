using Meshes

function create_centroids(gridsize::Int32, spfactor::Int32)
    # println("Build grid...")
    # grid properties
    lower_left = (0,0,0)
    N = Int64(gridsize)
    spcing_factor = Int64(spfactor)
    upper_right = (N,N,N) 
    spcing = (N*spcing_factor, N*spcing_factor, N*spcing_factor)

    # create 3D grid with N/N*spcing_factor spacing in each dimension, origin at (0,0,0)
    grid = Meshes.CartesianGrid(lower_left, upper_right, dims=spcing)
    # centroid of each cell
    # println("Build centroids...")
    centroids = Meshes.centroid.(grid)
    
    # reshape centroids 1D -> 3D 
    centroids = reshape([centroids...], N, N, N)

    return centroids
end