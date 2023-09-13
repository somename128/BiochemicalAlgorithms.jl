# using JLD2
# using Meshes, MeshViz
# using Makie, WGLMakie
using BenchmarkTools
using Profile
using ProfileView

include("correlation_docking.jl")
include("refine!.jl")
include("refine2!.jl")

res = Int32(1)
vdW = true
path_A = "testproteins/2hhb_alpha_chain.pdb"
path_B = "testproteins/2hhb_beta_chain.pdb"
correlation_docking(path_A, path_B, res, vdW)
@time sc = correlation_docking(path_A, path_B, res, vdW)

println(sc[1])
# save_object("src/dockings/results_docking.jld2", score)

# score_refined = refine2!(sc, Int32(10000), false)

# println(score_refined[1])
