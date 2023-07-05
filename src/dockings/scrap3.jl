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

#=
# do correlation_docking and store results 
# results are: scoring_table, grid_A, coord_B, centroids, gridsize
results_docking = correlation_docking("dummy_protein.pdb","dummy_ligand.pdb")

# save results in files
save_object("results_docking.jld2", results_docking)
=#
# load results of docking
results_docking = load_object("results_docking.jld2")

results_docking[1] == refine!(results_docking, Int32(100000))

