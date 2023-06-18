using BiochemicalAlgorithms
using FourierTools
using BenchmarkTools

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("extract_max.jl")
include("rotate_atoms.jl")

function generate_record(A::Array{Float32,3}, rotation::Matrix3{Float32}, roomcoordinates::Vector{Vector3{Float32}}, centroids::Array{Meshes.Point3, 3}, gridsize::Int32)
    # rotate atoms by rotation
    atoms = rotate_atoms(roomcoordinates, rotation, gridsize)
    # grid representation protein b
    B = grid_representation(atoms, gridsize, centroids)
    # fft-scoring
    C = ccorr(A,B; centered=false)
    # safe α,β,γ,R of max fft-scoring (c)
    max = extract_max(C)
    # build record for scoring table
    record = (α=max.α-1, β=max.β-1, γ=max.γ-1, R=rotation, score=max.score)

    return record
end