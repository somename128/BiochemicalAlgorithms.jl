using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

include("load_trans_pdb.jl")
include("helpers.jl")
include("min_max_atoms.jl")

sys = load_pdb("2ptc_protein.pdb")

molecules_df(sys).properties






