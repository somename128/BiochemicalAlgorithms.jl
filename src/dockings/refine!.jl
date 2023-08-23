using Rotations
using Distributions

include("generate_record.jl")
include("quaternion_functions.jl")
include("bingham_functions.jl")

function refine!(results_docking::Tuple, runs::Int32)

    # store scoring table, grid of protein A, roomcoordinates of protein b,
    # centroids and gridsize
    scoring_table = results_docking[1]
    grid_A = results_docking[2]
    roomcoordiantes_B = results_docking[3]
    centroids = results_docking[4]
    gridsize = results_docking[5]
    resolution = results_docking[6]

    # extract N best quaternions (current value: twenty)
    Q = Vector{QuaternionF32}()
    for i in eachrow(scoring_table)
        q = extract_quaternion(i)
        push!(Q, q)
    end

    # generate lookup table F for Bingham dirstribution
    F = create_lookup_table_F()
    # estimate V
    V = estimate_V(Q)
    # estimate Λ
    Λ = estimate_Λ(Q, V, F)
    # generate scatter matrix
    S = scatter_matrix(Q)

    # do refine runs number of times
    for i in 1:runs
        # sample new rotation
        R = metropolis_hastings_sampler(Q[1], Λ, V, S, F, Int32(100000))
        # generate record with new sampled rotation
        record = generate_record(grid_A, R, roomcoordiantes_B, centroids, gridsize, resolution)
        # check if record is better than first one in
        # current results
        if (results_docking[1][1,:].score[1] < record.score)
            push!(results_docking[1], record)
            # sort scoring_table
            sort!(results_docking[1], [:score], rev=[true])
        end
        # println(i,"/",runs)
    end

    # return new scoring table
    return results_docking[1]
end