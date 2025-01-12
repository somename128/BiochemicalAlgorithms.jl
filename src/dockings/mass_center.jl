function mass_center(atoms::Vector{Vector3{Float32}})
    # initalize vectors for storing x y z coordinates seperately
    X = Vector{Float32}()
    Y = Vector{Float32}()
    Z = Vector{Float32}()

    # fill vectors with coordinates
    for i in atoms
        push!(X, i[1])
        push!(Y, i[2])
        push!(Z, i[3])
    end

    # calculate center for every axis by sum up and divide by length
    mc_x = sum(X)/length(X)
    mc_y = sum(Y)/length(Y)
    mc_z = sum(Z)/length(Z)

    # store mass_center
    mass_center = Vector{Float32}([mc_x, mc_y, mc_z])

    return mass_center
end