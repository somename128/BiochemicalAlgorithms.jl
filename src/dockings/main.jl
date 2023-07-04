using BiochemicalAlgorithms
using BenchmarkTools
using ProfileView

include("correlation_docking.jl")

@time score = correlation_docking("src/dockings/dummy_protein_vol2.pdb","src/dockings/dummy_ligand.pdb")


