using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD2
using ProfileView
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations

include("correlation_docking.jl")
include("refine.jl")

# do correlation_docking
# score = correlation_docking("dummy_protein_vol2.pdb","dummy_ligand_vol3.pdb")
# save results
# save_object("score.jld2", score)

scoring_table = load_object("score.jld2")

refine(scoring_table,Int32(10))

