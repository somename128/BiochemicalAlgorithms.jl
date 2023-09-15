# using JLD2
# using Meshes, MeshViz
# using Makie, WGLMakie
using BenchmarkTools
using Profile
using ProfileView
using TimerOutputs
using JET

include("correlation_docking.jl")
include("refine!.jl")
include("refine2!.jl")

# for timetracking
# const to = TimerOutput()
res = Int32(1)
vdW = true
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
path_A = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
path_B = "src/dockings/testproteins/2hhb_beta_chain.pdb"
correlation_docking(init_A, init_B, res, vdW)
println("Initialization done.")
correlation_docking(path_A, path_B, res, vdW)
# show(to)
# println(sc[1])
# save_object("src/dockings/results_docking.jld2", score)

# score_refined = refine2!(sc, Int32(10000), false)

# println(score_refined[1])
