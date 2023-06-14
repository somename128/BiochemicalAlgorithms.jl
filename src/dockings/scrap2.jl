using BiochemicalAlgorithms

include("create_rotations.jl")
include("load_trans_pdb.jl")
include("extract_roomcoordinates.jl")
include("mass_center.jl")
include("rotate_atoms.jl")

N = Int32(64)
rotations = create_rotations()

protein_B = load_and_trans_pdb("src/dockings/dummy_ligand.pdb", N)
roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)

atoms = rotate_atoms(roomcoordiantes_atoms_B, rotations[1777], N)

length(rotations)
