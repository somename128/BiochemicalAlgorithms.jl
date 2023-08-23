# using JLD2
using Meshes, MeshViz
using Makie, WGLMakie

include("correlation_docking.jl")
include("refine!.jl")
include("refine2!.jl")

res = Int32(1)

@time score = correlation_docking("src/dockings/testproteins/2ptc_protein.pdb", "src/dockings/testproteins/2ptc_ligand.pdb", res, false)

println(score[1])
# save_object("src/dockings/results_docking.jld2", score)

@time score_refined = refine2!(refine2!(score, Int32(100)), Int32(100))

println(score_refined[1])
# viz(score_refined[2])
