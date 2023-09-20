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
using Meshes

#=
include("eval.jl")
include("eval_hhb.jl")

proteinA = "testproteins/2hhb_alpha_chain.pdb"
proteinB = "testproteins/2hhb_beta_chain.pdb"
complexAB = "testproteins/2hhb.pdb"

t = Vector3{Float32}(0, -1, -1)
R = (Float32(163.477), Float32(-13.0854), Float32(5.47113))

eva = eval_hhb(proteinA, proteinB, complexAB, R, t)
 
# CSV.write("equality.csv", DataFrame(eva), header = false)
=#

include("create_atomballs.jl")
include("load_trans_pdb.jl")
include("extract_roomcoordinates.jl")

protein_A = load_and_trans_pdb("testproteins/2hhb_alpha_chain.pdb", Int32(150))
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
atomballs = create_atomballs(roomcoordiantes_atoms_A, zero(Float32))
Int32(round(atomballs[1].center.coords[1]))