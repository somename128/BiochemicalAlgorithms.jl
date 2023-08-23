using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
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

include("evaluation.jl")

proteinA = "testproteins/2ptc_protein.pdb"
proteinB = "testproteins/2ptc_ligand.pdb"
complexAB = "testproteins/2ptc.pdb"
R = (Float32(-90), Float32(0), Float32(-90))
T = Vector3{Float32}(2, -1, 2)

evaluation(proteinA, proteinB, complexAB, R, T)
 