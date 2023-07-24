# using CSV

include("correlation_docking.jl")
include("refine!.jl")

res = Int32(2)

@time score = correlation_docking("src/dockings/dummy_protein.pdb", "src/dockings/dummy_ligand.pdb", res, true)

println(score[1])
# CSV.write("score_before_refinement.csv", score[1])

# @time score_refined = refine!(score, Int32(10))

# CSV.write("score_after_1000_refinement.csv", score_refined)