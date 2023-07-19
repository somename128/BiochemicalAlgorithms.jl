using DataFrames
using Base.Threads

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("create_rotations.jl")
include("create_centroids.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("set_gridsize.jl")

function correlation_docking(path_to_proteinA::String, path_to_proteinB::String) 
    # initialize scoring table
    scoring_table = DataFrame(α=zero(Float32), β=zero(Float32), γ=zero(Float32), R=(zero(Float32), zero(Float32), zero(Float32)), score=zero(Float32))
    # set gridsize N against protein size of greater protein
    N = set_gridsize(path_to_proteinA, path_to_proteinB)
    # load and translate protein a
    protein_A = load_and_trans_pdb(path_to_proteinA, N)
    roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
    # load and translate protein b
    protein_B = load_and_trans_pdb(path_to_proteinB, N)
    roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
    # calculate centroids in a NxNxN grid with cells
    # of 1 angström, only done once
    centroids = create_centroids(N, one(Int32))
    # grid representation protein a
    A = grid_representation(roomcoordiantes_atoms_A, N, centroids, false)
    # get rotations - build via rigid_transform! with translation vector (0,0,0)
    rotations = create_rotations()
    # lock for threads (unsure how this really works)
    lk = ReentrantLock() 
    # rotate protein b by R
    @threads for i in eachindex(rotations)
        # generate record for scoring table
        record = generate_record(A, rotations[i], roomcoordiantes_atoms_B, centroids,N)
        lock(lk) do
            #=
            if(scoring_table.score[1] < record.score)
                scoring_table[1,:] = (α=record.α, β=record.β, γ=record.γ, R=record.R, score=record.score)
                println(i)
            end
            =#
            push!(scoring_table, record)
        end
        # println(i,"/",length(rotations))
    end

    # sort scoring_table
    sort!(scoring_table, [:score], rev=[true])

    # return five greatest values, grid representation of A,
    # roomcoordinates of B, centroids and gridsize for refinement
    return scoring_table[1:20, :], A, roomcoordiantes_atoms_B, centroids, N
end
