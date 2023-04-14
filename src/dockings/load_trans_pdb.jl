using BiochemicalAlgorithms

function load_and_trans_pdb(path_to_pdb)

    # load protein data from PDB
    println("Load PDB file...")
    protein = load_pdb(path_to_pdb)

    # translate protein in positive space
    println("Translate protein in positive space...")
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

    # calculate x y z min for translation (max maybe needed one time)
    min_x = minimum(X)
    # max_x = maximum(X)
    min_y = minimum(Y)
    # max_y = maximum(Y)
    min_z = minimum(Z)
    # max_z = maximum(Z)

    # set translation vector depending on atom coordinates
    if (min_x < 0 && min_y < 0 && min_z < 0)
        translation_vector = Vector3{Float32}(abs(min_x)+1, abs(min_y)+1, abs(min_z)+1)
    elseif (min_x < 0 && min_y < 0 && min_z >= 0)
        translation_vector = Vector3{Float32}(abs(min_x)+1, abs(min_y)+1, 0)
    elseif (min_x < 0 && min_y >= 0 && min_z < 0)
        translation_vector = Vector3{Float32}(abs(min_x)+1, 0, abs(min_z)+1)   
    elseif (min_x < 0 && min_y >= 0 && min_z >= 0)
        translation_vector = Vector3{Float32}(abs(min_x)+1, 0, 0)
    elseif (min_x >= 0 && min_y < 0 && min_z < 0)
        translation_vector = Vector3{Float32}(0, abs(min_y)+1, abs(min_z)+1)
    elseif (min_x >= 0 && min_y < 0 && min_z >= 0)
        translation_vector = Vector3{Float32}(0, abs(min_y)+1, 0) 
    elseif (min_x >= 0 && min_y >= 0 && min_z < 0)
        translation_vector = Vector3{Float32}(0, 0, abs(min_z)+1)
    else 
        translation_vector = (0,0,0)
    end

    # translate protein in "positive space"
    BiochemicalAlgorithms.translate!(protein,translation_vector)

    return protein
end
