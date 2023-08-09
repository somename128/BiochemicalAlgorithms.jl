using BiochemicalAlgorithms
using Meshes
using Quaternions

include("load_trans_pdb.jl")
include("extract_roomcoordinates.jl")
include("rotate_atoms.jl")
include("create_rotations.jl")

N = Int32(128)
protein = load_and_trans_pdb("dummy_protein_vol4.pdb", N)
atoms = extract_roomcoordinates(protein)
# radii = Dict("C" => 1.7, "H" => 1.0, "N" => 1.5, "O" => 1.4)
# rotations = create_rotations() 
# rotated_tuple = rotate_atoms(atoms, rotations[2], N)
