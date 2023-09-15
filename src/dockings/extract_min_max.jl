function extract_min_max(roomcoordinates::Vector{Vector3{Float32}})

    # initalize vectors for storing x y z coordinates seperately
    X = Array{Float32}(undef, length(roomcoordinates))
    Y = Array{Float32}(undef, length(roomcoordinates))
    Z = Array{Float32}(undef, length(roomcoordinates))

    # fill vectors with coordinates
    for i in eachindex(roomcoordinates)
        X[i] = roomcoordinates[i][1]
        Y[i] = roomcoordinates[i][2]
        Z[i] = roomcoordinates[i][3]
    end

    # calculate min max per axis 
    min_x = minimum(X)
    max_x = maximum(X)
    min_y = minimum(Y)
    max_y = maximum(Y)
    min_z = minimum(Z)
    max_z = maximum(Z)

    return min_x, max_x, min_y, max_y, min_z, max_z
end