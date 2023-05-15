using BiochemicalAlgorithms
using BenchmarkTools

include("correlation_docking.jl")

score = correlation_docking("2ptc_protein.pdb","2ptc_ligand.pdb", 128)

score[findmax(score.score)[2]]
