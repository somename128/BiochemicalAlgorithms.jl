# using JLD2
using Meshes, MeshViz
using Makie, WGLMakie

include("correlation_docking.jl")
include("refine!.jl")
include("refine2!.jl")

res = Int32(1)

@time score = correlation_docking("src/dockings/dummy_protein.pdb", "src/dockings/dummy_ligand.pdb", res, false)

println(score[1])
# save_object("src/dockings/results_docking.jld2", score)

@time score_refined = refine!(score, Int32(10000))

println(score_refined)
# viz(score_refined[2])
