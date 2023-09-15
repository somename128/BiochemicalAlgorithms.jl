function mass_center(atoms::Vector{Vector3{Float32}})
    # initalize vectors for storing x y z coordinates seperately
    X = Array{Float32}(undef, length(atoms))
    Y = Array{Float32}(undef, length(atoms))
    Z = Array{Float32}(undef, length(atoms))

    # fill vectors with coordinates
    for i in eachindex(atoms)
        X[i] =  atoms[i][1]
        Y[i] =  atoms[i][2]
        Z[i] =  atoms[i][3]
    end

    # calculate center for every axis by sum up and divide by length
    mc_x = sum(X)/length(X)
    mc_y = sum(Y)/length(Y)
    mc_z = sum(Z)/length(Z)

    # store mass_center
    mass_center = Vector{Float32}([mc_x, mc_y, mc_z])

    return mass_center
end