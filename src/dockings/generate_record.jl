using BiochemicalAlgorithms
using FFTW
using BenchmarkTools

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("extract_max.jl")

function generate_record(A::Array{Float64,3}, rotation::RigidTransform{Float32}, path_to_proteinB::String, centroids::Array{Meshes.Point3, 3}, gridsize::Int64)
    # load ant trans pdb each time to get origin coordinates
    protein = load_and_trans_pdb(path_to_proteinB, gridsize)
    # rotate protein by rotation
    rigid_transform!(protein, rotation)
    # grid representation protein b
    B = grid_representation(protein, gridsize, centroids)
    # fft-scoring
    C = ifft(fft(A).*fft(B))
    # safe α,β,γ,R of max fft-scoring (c)
    max = extract_max(C)
    # build record for scoring table
    record = (α=max.α, β=max.β, γ=max.γ, R=rotation, score=max.score)

    return record
end