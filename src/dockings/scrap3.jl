using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

include("load_trans_pdb.jl")
include("helpers.jl")
include("min_max_atoms.jl")
include("mass_center.jl")
include("create_centroids.jl")
include("extract_roomcoordinates.jl")
include("create_atomballs.jl")
include("create_rotations.jl")
include("mass_center.jl")
include("rotate_atoms.jl")

N = 128
protein = load_and_trans_pdb("2ptc_ligand.pdb", N)
# atomballs = create_atomballs(protein)
rotations = create_rotations()
# @time mc = mass_center(protein)
# @time min_max = min_max_atoms(protein)

# atoms_m = Vector{Vector3{Float32}}()
atoms = extract_roomcoordinates(protein)

f(v::Vector3{Float32}) = rotations[5] * v 

atoms_f = f.(atoms)

atoms_r = rotate_atoms(atoms, rotations[5])

println(atoms ≈ atoms_f)
println(atoms_r ≈ atoms_f)





