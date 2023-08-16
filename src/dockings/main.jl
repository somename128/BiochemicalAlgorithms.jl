# using JLD2
using Meshes, MeshViz
using Makie, WGLMakie

include("correlation_docking.jl")
include("refine!.jl")

res = Int32(1)

@time score = correlation_docking("src/dockings/dummy_protein_vol4.pdb", "src/dockings/dummy_ligand_vol4.pdb", res, true)

println(score[1])
# save_object("src/dockings/results_docking.jld2", score)

@time score_refined = refine!(score, Int32(1000))

println(score_refined[1])
# viz(score_refined[2])
