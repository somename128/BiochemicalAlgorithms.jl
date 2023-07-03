using BiochemicalAlgorithms
using BenchmarkTools
using ProfileView

include("correlation_docking.jl")

@time score = correlation_docking("src/dockings/dummy_protein.pdb","src/dockings/dummy_ligand_vol2.pdb")


