using JLD2
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

res = Int32(1)
# grid surface, 120-cell
vdW = false
hypercube = true
# 2hhb
println("2hhb")
init = true
init_N = Int32(150)
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
correlation_docking(init_A, init_B, res, vdW, hypercube, init, init_N)
println("Initialization done.")
init = false
path_A = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
path_B = "src/dockings/testproteins/2hhb_beta_chain.pdb"
@time score = correlation_docking(path_A, path_B, res, vdW, hypercube, init, init_N)
save_object("src/dockings/testrun/2hhb_grid_120.jld2", score[1])
save_object("src/dockings/testrun_huge/2hhb_grid_120.jld2", score)
println("\n")
# 2mhb
println("2mhb")
init = true
init_N = Int32(128)
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
correlation_docking(init_A, init_B, res, vdW, hypercube, init, init_N)
println("Initialization done.")
init = false
path_A = "src/dockings/testproteins/2mhb_alpha_chain.pdb"
path_B = "src/dockings/testproteins/2mhb_beta_chain.pdb"
@time score = correlation_docking(path_A, path_B, res, vdW, hypercube, init, init_N)
save_object("src/dockings/testrun/2mhb_grid_120.jld2", score[1])
save_object("src/dockings/testrun_huge/2mhb_grid_120.jld2", score)
println("\n")
# 2ptc
println("2ptc")
init = true
init_N = Int32(128)
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
correlation_docking(init_A, init_B, res, vdW, hypercube, init, init_N)
println("Initialization done.")
init = false
path_A = "src/dockings/testproteins/2ptc_proteinase.pdb"
path_B = "src/dockings/testproteins/2ptc_inhibitor.pdb"
@time score = correlation_docking(path_A, path_B, res, vdW, hypercube, init, init_N)
save_object("src/dockings/testrun/2ptc_grid_120.jld2", score[1])
save_object("src/dockings/testrun_huge/2ptc_grid_120.jld2", score)
println("\n")
# 3apr
println("3apr")
init = true
init_N = Int32(128)
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
correlation_docking(init_A, init_B, res, vdW, hypercube, init, init_N)
println("Initialization done.")
init = false
path_A = "src/dockings/testproteins/3apr_proteinase.pdb"
path_B = "src/dockings/testproteins/3apr_inhibitor.pdb"
@time score = correlation_docking(path_A, path_B, res, vdW, hypercube, init, init_N)
save_object("src/dockings/testrun/3apr_grid_120.jld2", score[1])
save_object("src/dockings/testrun_huge/3apr_grid_120.jld2", score)
println("\n")
# 3ts1
println("3ts1")
init = true
init_N = Int32(150)
init_A = "src/dockings/simple_geometry/cube_origin.pdb"
init_B = "src/dockings/simple_geometry/cube_origin.pdb"
correlation_docking(init_A, init_B, res, vdW, hypercube, init, init_N)
println("Initialization done.")
init = false
path_A = "src/dockings/testproteins/3ts1_protein.pdb"
path_B = "src/dockings/testproteins/3ts1_ligand.pdb"
@time score = correlation_docking(path_A, path_B, res, vdW, hypercube, init, init_N)
save_object("src/dockings/testrun/3ts1_grid_120.jld2", score[1])
save_object("src/dockings/testrun_huge/3ts1_grid_120.jld2", score)
println("\n")
