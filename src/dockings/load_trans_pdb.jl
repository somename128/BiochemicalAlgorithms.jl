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

    # calculate x y z min and max for translation that 
    # needed for optimal rotation
    min_x = minimum(X)
    max_x = maximum(X)
    min_y = minimum(Y)
    max_y = maximum(Y)
    min_z = minimum(Z)
    max_z = maximum(Z)

    # calculate "largest" diameter of protein
    # translation adding this value in x y z ensure that
    # protein do not "vanish" through rotation
    # so far: for maximum function vector needed
    RAD = Vector{Float32}()
    push!(DIA,abs(min_x-max_x)/2)
    push!(DIA,abs(min_y-max_y)/2)
    push!(DIA,abs(min_z-max_z)/2)

    # maximum radius of x y z added 3 for atom radians purposes
    max_rad = maximum(DIA)+3

    # set translation vector depending on atom coordinates
    # and maximum radius in x y z 
    if (min_x < 0 && min_y < 0 && min_z < 0)
        translation_vector = Vector3{Float32}(abs(min_x)+max_rad, 
            abs(min_y)+max_rad, abs(min_z)+max_rad)
    elseif (min_x < 0 && min_y < 0 && min_z >= 0)
        translation_vector = Vector3{Float32}(abs(min_x)+max_rad, 
            abs(min_y)+max_rad, 0+max_rad)
    elseif (min_x < 0 && min_y >= 0 && min_z < 0)
        translation_vector = Vector3{Float32}(abs(min_x)+max_rad, 
            0+max_rad, abs(min_z)+max_rad)   
    elseif (min_x < 0 && min_y >= 0 && min_z >= 0)
        translation_vector = Vector3{Float32}(abs(min_x)+max_rad, 
            0+max_rad, 0+max_rad)
    elseif (min_x >= 0 && min_y < 0 && min_z < 0)
        translation_vector = Vector3{Float32}(0+max_rad, 
            abs(min_y)+max_rad, abs(min_z)+max_rad)
    elseif (min_x >= 0 && min_y < 0 && min_z >= 0)
        translation_vector = Vector3{Float32}(0+max_rad, 
            abs(min_y)+max_rad, 0+max_rad) 
    elseif (min_x >= 0 && min_y >= 0 && min_z < 0)
        translation_vector = Vector3{Float32}(0+max_rad, 
            0+max_rad, abs(min_z)+max_rad)
    else 
        translation_vector = (0+max_rad,0+max_rad,0+max_rad)
    end

    # translate protein in "positive space"
    BiochemicalAlgorithms.translate!(protein,translation_vector)

    return protein
end
