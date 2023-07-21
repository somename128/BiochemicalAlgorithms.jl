using Distributions
using Rotations
using Quaternions

include("set_marked_cells.jl")
include("load_trans_pdb.jl")
include("create_atomballs.jl")
include("create_centroids.jl")
include("min_max_atoms.jl")
include("extract_roomcoordinates.jl")

N = Int32(64)
resolution = Int32(2)
protein_A = load_and_trans_pdb("dummy_protein.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
centroids = create_centroids(N, res)
atomballs = create_atomballs(roomcoordiantes_atoms_A, Float32(1))
# marked_cells = set_marked_cells(atomballs, centroids, roomcoordiantes_atoms_A, res)
min_max = min_max_atoms(roomcoordiantes_atoms_A)
centroids[1]