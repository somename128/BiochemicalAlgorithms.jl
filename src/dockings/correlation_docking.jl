using DataFrames
using Base.Threads

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("create_rotations.jl")
include("create_centroids.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("set_gridsize.jl")
include("create_rotations2.jl")

# function that takes two paths to pdb structure of proteins and do correlation docking
# it also gets a boolean if vdW surface is choosen
# additionally it gets the resolution which means how small the cubes of the
# grid representation should be
function correlation_docking(path_to_proteinA::String, path_to_proteinB::String, resolution::Int32, vdW::Bool)
    # initialize scoring table
    scoring_table = DataFrame(α=zero(Float32), β=zero(Float32), γ=zero(Float32), R=(zero(Float32), zero(Float32), zero(Float32)), score=zero(Float32))
    # set gridsize N against protein size of greater protein
    # N = set_gridsize(path_to_proteinA, path_to_proteinB)
    N = Int32(32)
    # load and translate protein a
    protein_A = load_and_trans_pdb(path_to_proteinA, N)
    roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
    # load and translate protein b
    protein_B = load_and_trans_pdb(path_to_proteinB, N)
    roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
    # calculate centroids in a NxNxN grid with cells
    # of 1 angström, only done once
    centroids = create_centroids(N, resolution)
    # grid representation protein a
    A = grid_representation(roomcoordiantes_atoms_A, N, centroids, resolution, false, vdW)
    # get quaternion rotations (via 20 degree or 120-cell)
    rotations = create_rotations2()
    # lock for threads (unsure how this really works)
    lk = ReentrantLock() 
    # rotate protein b by R
    @threads for i in eachindex(rotations)
        # generate record for scoring table
        record = generate_record(A, rotations[i], roomcoordiantes_atoms_B, centroids, N, resolution, vdW)
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
    return scoring_table[1:10, :], A, roomcoordiantes_atoms_B, centroids, N, resolution
end
