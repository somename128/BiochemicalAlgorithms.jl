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

include("evaluation_simple.jl")

proteinA = "simple_geometry/cube_origin_huge_A.pdb"
proteinB = "simple_geometry/cube_origin_huge_B.pdb"
complexAB = "simple_geometry/cube_origin_huge.pdb"

R = (Float32(0), Float32(0), Float32(0))
T = Vector3{Float32}(-2, -2, -2)

evaluation_simple(proteinA, proteinB, complexAB, R, T)
 
# CSV.write("equality.csv", DataFrame(eva), header = false)