using BiochemicalAlgorithms

include("mass_center.jl")
include("extract_roomcoordinates.jl")

function load_and_trans_pdb(path_to_pdb::String, gridsize::Int32)

    # load protein data from PDB
    # first system, then molecule
    # datatype protein: Molecule{Float32}
    # println("Load PDB file...")
    protein = molecules(load_pdb(path_to_pdb))[1]
    atoms = [i.r for i in eachrow(atoms_df(protein))]
    # translate protein in center of grid
    # println("Translate protein in center of grid...")
    # get mass_center of protein
    mc = mass_center(atoms)

    # set translation vector depending on mass center
    # set mass center in center of grid
    # multiplication by -1 to change sign,
    # translation vector is the "inverse" of mass_center minus
    # grid_center (in this case N = 128, so (64,64,64)) 
    translation_vector = (-1)*Vector3{Float32}(mc[1]-gridsize/2, mc[2]-gridsize/2, 
        mc[3]-gridsize/2)

    # translate protein in to center of grid
    BiochemicalAlgorithms.translate!(protein,translation_vector)

    return protein
end
