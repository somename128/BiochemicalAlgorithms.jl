using BiochemicalAlgorithms
using BenchmarkTools
using JLD
using TypedTables

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("create_rotations.jl")
include("create_centroids.jl")
include("generate_record.jl")

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
    # get rotations - build via rigid_transform! with translation vector (0,0,0)
    rotations = create_rotations()

    # rotate protein b by R
    for i in eachindex(rotations)
        # generate record for scoring table
        record = generate_record(A,rotations[i],path_to_proteinB,centroids,N)
        push!(scoring_table,record)
        println(i)
    end
    # --------------------------------------------
    # safe results of run
    # save("scoring_table.jld", "scoring_table", scoring_table)
    # --------------------------------------------

    # loop to rotation if "greatest" c not reached
    return scoring_table
end
