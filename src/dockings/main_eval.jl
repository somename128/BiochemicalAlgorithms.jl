using JLD2
using BiochemicalAlgorithms
using Meshes, MeshViz
using Makie, WGLMakie
using BenchmarkTools
using Profile
using ProfileView
# using TimerOutputs
# using JET
using WAV
using Plots

# include("correlation_docking.jl")
# include("refine!.jl")
# include("refine2!.jl")
include("evalscript.jl")
include("eval_hhb.jl")
include("eval.jl")

# load protein A,B and complex AB
# --------------
# |CHANGE HERE!|
# --------------
pathA = "src/dockings/testproteins/3ts1_protein.pdb"
pathB = "src/dockings/testproteins/3ts1_ligand.pdb"
pathC = "src/dockings/testproteins/3ts1.pdb"
proteinA = molecules(load_pdb(pathA))[1]
proteinB = molecules(load_pdb(pathB))[1]
complexAB = molecules(load_pdb(pathC))[1]

# load data, insert new column
# --------------
# |CHANGE HERE!|
# --------------
sc = load_object("src/dockings/testrun_huge/3ts1_grid_120.jld2")
insertcols!(sc[1], :rmsd => Float32[typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32)])

for i in eachrow(sc[1])
    # calculate rmsd for record
    local t = Vector3{Float32}(i.α, i.β, i.γ)
    local R = (i.R[1], i.R[2], i.R[3])
    # --------------
    # |CHANGE HERE!|
    # --------------
    i.rmsd = try eval(proteinA, proteinB, complexAB, R, t)
    catch
        typemax(Float32)
    end
end
# show(stdout, MIME("text/latex"), sc[1])

# --------------
# |CHANGE HERE!|
# --------------
vdW = false
λ = Float32(100)
# plot refinement timeseries
# --------------
# |CHANGE HERE!|
# --------------
rmsds = evalscript(proteinA, proteinB, complexAB, deepcopy(sc), λ, vdW)
Plots.plot(rmsds, xlabel="iteration", ylabel="rmsd", ylims=(0, 25), xticks=0:10:length(rmsds), label="1", title="Refinement 3TS1, vdW, 120-cell", linewidth=3)
rmsds = evalscript(proteinA, proteinB, complexAB, deepcopy(sc), λ, vdW)
Plots.plot!(rmsds, xlabel="iteration", ylabel="rmsd", ylims=(0, 25), xticks=0:10:length(rmsds), label="2", title="Refinement 3TS1, vdW, 120-cell", linewidth=3)
rmsds = evalscript(proteinA, proteinB, complexAB, deepcopy(sc), λ, vdW)
Plots.plot!(rmsds, xlabel="iteration", ylabel="rmsd", ylims=(0, 25), xticks=0:10:length(rmsds), label="3", title="Refinement 3TS1, vdW, 120-cell", linewidth=3)
rmsds = evalscript(proteinA, proteinB, complexAB, deepcopy(sc), λ, vdW)
Plots.plot!(rmsds, xlabel="iteration", ylabel="rmsd", ylims=(0, 25), xticks=0:10:length(rmsds), label="4", title="Refinement 3TS1, vdW, 120-cell", linewidth=3)
rmsds = evalscript(proteinA, proteinB, complexAB, deepcopy(sc), λ, vdW)
Plots.plot!(rmsds, xlabel="iteration", ylabel="rmsd", ylims=(0, 25), xticks=0:10:length(rmsds), label="5", title="Refinement 3TS1, vdW, 120-cell", linewidth=3)

# save plot
# --------------
# |CHANGE HERE!|
# --------------
savefig("src/dockings/latex/3ts1/refine_runs_vdW_120.png") 

y, fs = wavread(raw"src/dockings/ff_victory.wav")
wavplay(y, fs)