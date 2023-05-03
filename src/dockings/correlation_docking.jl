using BiochemicalAlgorithms
using BenchmarkTools
using FFTW
using JLD
using TypedTables

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("extract_max.jl")
include("create_rotations.jl")
include("create_centroids.jl")

function correlation_docking(path_to_proteinA::String, path_to_proteinB::String, gridsize::Int64)
    # initialize scoring table
    scoring_table = Table(α=[], β=[], γ=[], R=[], score=[])
    # set grid size N 
    N = gridsize
    # load and translate protein a
    protein_A = load_and_trans_pdb(path_to_proteinA, N)
    # load and translate protein b
    # protein_B = load_and_trans_pdb(path_to_proteinB)
    # calculate centroids in a NxNxN grid with cells
    # of 1 angström, only done once
    centroids = create_centroids(N, 1)
    # grid representation protein a
    A = grid_representation(protein_A, N, centroids)
    # save("grid_repo_A.jld", "A", A)
    # A = load("grid_rep_A.jld")["A"]

    # get rotations - build via rigid_transform! with translation vector (0,0,0)
    rotations = create_rotations()
    # rotate protein b by R
    @time begin 
    for i in eachindex(rotations)
        # load and translate protein b per loop to use rotations relative to 
        # origin coordinates
        protein_B = load_and_trans_pdb(path_to_proteinB, N)
        rigid_transform!(protein_B, rotations[i])
        # grid representation protein b
        B = grid_representation(protein_B, N, centroids)
        # fft-scoring
        C = ifft(fft(A).*fft(B))
        # -----------------------------------------
        # only for time saving purposes, not part of the real algorithm
        # save("C_matrix.jld", "C", C)

        # C = load("C_matrix.jld")["C"]
        # -----------------------------------------
        # safe α,β,γ,R of max fft-scoring (c)
        max = extract_max(C)
        # build record for scoring table
        record = (α=max.α, β=max.β, γ=max.γ, R=rotations[i], score=max.score)
        push!(scoring_table,record)
    end
    end    
    # --------------------------------------------
    # safe results of run
    # save("scoring_table.jld", "scoring_table", scoring_table)
    # --------------------------------------------

    # loop to rotation if "greatest" c not reached
    return scoring_table
end
