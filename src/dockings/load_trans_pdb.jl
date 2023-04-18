using BiochemicalAlgorithms

include("mass_center.jl")

function load_and_trans_pdb(path_to_pdb)

    # load protein data from PDB
    println("Load PDB file...")
    protein = load_pdb(path_to_pdb)

    # translate protein in center of grid
    println("Translate protein in center of grid...")
    # get mass_center of protein
    mc = mass_center(protein)

    # set translation vector depending on mass center
    # set mass center in center of grid
    # multiplication by -1 to change sign,
    # translation vector is the "inverse" of mass_center minus
    # grid_center (in this case N = 128, so (64,64,64)) 
    translation_vector = (-1)*Vector3{Float32}(mc[1]-64, mc[2]-64, 
        mc[3]-64)

    # translate protein in to center of grid
    BiochemicalAlgorithms.translate!(protein,translation_vector)

    return protein
end
