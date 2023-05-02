using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

include("load_trans_pdb.jl")
include("helpers.jl")
include("min_max_atoms.jl")

protein = load_and_trans_pdb("2ptc_ligand.pdb") 

min_max_atoms(protein)[2]





