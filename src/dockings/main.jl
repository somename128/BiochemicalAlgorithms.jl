# using CSV

include("correlation_docking.jl")
include("refine!.jl")

@time score = correlation_docking("src/dockings/dummy_protein.pdb", "src/dockings/dummy_ligand.pdb")

println(score[1])
# CSV.write("score_before_refinement.csv", score[1])

# @time score_refined = refine!(score, Int32(1000))

# CSV.write("score_after_1000_refinement.csv", score_refined)