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
