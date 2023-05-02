using BiochemicalAlgorithms

include("correlation_docking.jl")

@time begin 
score = correlation_docking("2ptc_protein.pdb","2ptc_ligand.pdb", 128)

score[findmax(score.score)[2]]
end