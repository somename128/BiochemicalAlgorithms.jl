using BiochemicalAlgorithms
using BenchmarkTools
using JLD2
using ProfileView
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations
using DataFrames
using DelimitedFiles
using CSV
using Meshes, MeshViz
using Makie, WGLMakie
using Plots
using WAV
using PDBTools

include("eval_hhb.jl")
include("eval.jl")
include("quaternion_functions.jl")

#=
pathA = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
proteinA = molecules(load_pdb(pathA))[1]
for i in eachrow(atoms_df(proteinA))
    println(i.r[1])
end
=#

protein = readPDB("src/dockings/ballview/2hhb.pdb")

for i in eachrow(protein)
    protein.i[1].x
end

# sound
#=
y, fs = wavread(raw"src/dockings/ff_victory.wav")
wavplay(y, fs)
=#