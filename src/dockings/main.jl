using BiochemicalAlgorithms
using BenchmarkTools
using ProfileView

include("correlation_docking.jl")

@profview score = correlation_docking("2ptc_protein.pdb","2ptc_ligand.pdb", 128)


