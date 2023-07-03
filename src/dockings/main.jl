using BiochemicalAlgorithms
using BenchmarkTools
using ProfileView

include("correlation_docking.jl")

@time score = correlation_docking("src/dockings/2ptc_protein.pdb","src/dockings/2ptc_ligand.pdb")


