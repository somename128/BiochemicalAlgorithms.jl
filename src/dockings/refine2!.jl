using Rotations
using Distributions
using Base.Threads

include("generate_record.jl")
include("quaternion_functions.jl")
include("bingham_functions_vol2.jl")

function refine2!(results_docking::Tuple, runs::Int32, vdW::Bool)

    # store scoring table, grid of protein A, roomcoordinates of protein b,
    # centroids and gridsize
    scoring_table = results_docking[1]
    grid_A = results_docking[2]
    roomcoordiantes_B = results_docking[3]
    centroids = results_docking[4]
    gridsize = results_docking[5]
    resolution = results_docking[6]

    # extract N best quaternions (current value: twenty)
    # TODO: only get best rotation
    Q = Vector{QuaternionF32}()
    for i in eachrow(scoring_table)
        q = extract_quaternion(i)
        push!(Q, q)
    end

    # set parameters for quaternion sampling
    μ = [Q[1].s, Q[1].v1, Q[1].v2, Q[1].v3]
    λ = 2

    # sample quaternions Bingham distributed around best quaternion
    rotations = sample_quaternions(μ, λ, runs)
    # lock for threads (unsure how this really works)
    lk = ReentrantLock()
    # calculate scorings for sampled rotations
    @threads for i in eachindex(rotations)
        # generate rotation quaternion
        R = QuaternionF32(rotations[i][1], rotations[i][2], rotations[i][3], rotations[i][4])
        # generate record with new sampled rotation
        record = generate_record(grid_A, R, roomcoordiantes_B, centroids, gridsize, resolution, vdW)
        # check if record is better than first one in
        # current results
        lock(lk) do
            if (results_docking[1][1,:].score[1] < record.score)
                push!(results_docking[1], record)
                # sort scoring_table
                sort!(results_docking[1], [:score], rev=[true])
            end
        end
        # println(i,"/",runs)
    end

    # return new scoring table
    return results_docking
end