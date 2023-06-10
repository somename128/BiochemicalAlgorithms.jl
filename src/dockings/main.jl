using BiochemicalAlgorithms
using BenchmarkTools
using ProfileView

include("correlation_docking.jl")

@time score = correlation_docking("src/dockings/dummy_protein.pdb","src/dockings/dummy_ligand.pdb", Int32(128))


