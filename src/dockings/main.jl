using JLD2
# using Meshes, MeshViz
# using Makie, WGLMakie
using BenchmarkTools
using Profile
using ProfileView
using TimerOutputs
using JET
using WAV

include("correlation_docking.jl")
include("refine!.jl")
include("refine2!.jl")

#=
# for timetracking
# const to = TimerOutput()
res = Int32(1)
vdW = false
hypercube = true
init = false
init_N = Int(32)
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
path_A = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
path_B = "src/dockings/testproteins/2hhb_beta_chain.pdb"
correlation_docking(init_A, init_B, res, vdW, hypercube, init, init_N)
println("Initialization done.")
@time correlation_docking(path_A, path_B, res, vdW, hypercube, init, init_N)
# show(to)
# println(sc[1])
# save_object("src/dockings/results_docking.jld2", score)
=# 

sc = load_object("src/dockings/testrun_huge/2hhb_vdW_120.jld2")
println(sc[1])
# refine2!(sc, Int32(10), false)
@time score_refined = refine2!(sc, Int32(50), true)
println(score_refined[1][1:10, :])

@time score_refined = refine2!(sc, Int32(100), true)
println(score_refined[1][1:10, :])

@time score_refined = refine2!(sc, Int32(150), true)
println(score_refined[1][1:10, :])

y, fs = wavread(raw"src/dockings/ff_victory.wav")
wavplay(y, fs)
