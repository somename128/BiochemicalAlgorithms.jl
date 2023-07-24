using Distributions
using Rotations
using Quaternions

include("set_marked_cells.jl")
include("load_trans_pdb.jl")
include("create_atomballs.jl")
include("create_centroids.jl")
include("min_max_atoms.jl")
include("extract_roomcoordinates.jl")
include("quaternion_functions.jl")
include("create_rotations.jl")

#=
N = Int32(64)
resolution = Int32(2)
protein_A = load_and_trans_pdb("dummy_protein.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
centroids = create_centroids(N, res)
atomballs = create_atomballs(roomcoordiantes_atoms_A, Float32(1))
# marked_cells = set_marked_cells(atomballs, centroids, roomcoordiantes_atoms_A, res)
min_max = min_max_atoms(roomcoordiantes_atoms_A)
=#
rotations = create_rotations()
vs = Vector{Vector3{Float32}}()
a = Vector3{Float32}(1,2,3)
b = Vector3{Float32}(4,5,6)
push!(vs,a)
push!(vs,b)
g(v::Vector3{Float32}) = rotate_vector(rotations[11], v)
b = g.(vs)

