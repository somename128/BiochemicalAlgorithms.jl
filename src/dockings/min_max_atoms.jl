using BiochemicalAlgorithms

function min_max_atoms(protein)
    # extract room coordinates of atoms of the protein
    atoms_in_space = protein.atoms.r
    # initalize vectors for storing x y z coordinates seperately
    X = Vector{Float32}()
    Y = Vector{Float32}()
    Z = Vector{Float32}()

    # fill vectors with coordinates
    for i in atoms_in_space
        push!(X, i[1])
        push!(Y, i[2])
        push!(Z, i[3])
    end

    # calculate min max per axis and round+int it (for relevant
    # centroids indexing)
    min_x = floor(Int,minimum(X))
    max_x = ceil(Int,maximum(X))
    min_y = floor(Int,minimum(Y))
    max_y = ceil(Int,maximum(Y))
    min_z = floor(Int,minimum(Z))
    max_z = ceil(Int,maximum(Z))

    return min_x, max_x, min_y, max_y, min_z, max_z
end