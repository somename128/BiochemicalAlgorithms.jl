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

include("bingham_functions_vol2.jl")
include("vmf.jl")

q = Float32[0.235, 0.568, 0.323, 0.6845]

# sample_quaternions(q, Int32(2), Int32(10))
for i in 1:10
    println(sampleVMF(q, 2))
end