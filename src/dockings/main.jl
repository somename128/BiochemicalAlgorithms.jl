# using JLD2
# using Meshes, MeshViz
# using Makie, WGLMakie
using ProfileView

include("correlation_docking.jl")
include("refine!.jl")
include("refine2!.jl")

res = Int32(2)

score = correlation_docking("src/dockings/simple_geometry/dot_origin_vdW.pdb", "src/dockings/simple_geometry/dot_origin_vdW.pdb", res, true)

println(score[1])
# save_object("src/dockings/results_docking.jld2", score)

score_refined = refine2!(score, Int32(10000), true)

println(score_refined[1])
# viz(score_refined[2])
