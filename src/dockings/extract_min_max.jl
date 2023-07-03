function extract_min_max(roomcoordinates::Vector{Vector3{Float32}})

    # initalize vectors for storing x y z coordinates seperately
    X = Vector{Float32}()
    Y = Vector{Float32}()
    Z = Vector{Float32}()

    # fill vectors with coordinates
    for i in roomcoordinates
        push!(X, i[1])
        push!(Y, i[2])
        push!(Z, i[3])
    end

    # calculate min max per axis and round+int it (for relevant
    # centroids indexing)
    # adding -2 and +2 to add extra space to get all relevant
    # centroids which are inside atom radius (=1.8A) 
    # if a value is less than one it's set to one
    min_x = minimum(X)
    max_x = maximum(X)
    min_y = minimum(Y)
    max_y = maximum(Y)
    min_z = minimum(Z)
    max_z = maximum(Z)

    return min_x, max_x, min_y, max_y, min_z, max_z



end