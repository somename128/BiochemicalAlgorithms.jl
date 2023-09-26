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


include("bingham_functions_vol2.jl")

μ = Float32[0.234, 0.345, 0.235, -0.46347]
λ = Float32(100)

sample_quaternions(μ, λ, Int32(10))
