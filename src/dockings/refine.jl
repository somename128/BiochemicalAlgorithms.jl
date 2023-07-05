using Rotations
using Distributions

include("generate_record.jl")

function refine!(results_docking::Tuple, runs::Int32)

    # store grid of protein A, roomcoordinates of protein b,
    # centroids and gridsize
    grid_A = results_docking[2]
    roomcoordiantes_B = results_docking[3]
    centroids = results_docking[4]
    gridsize = results_docking[5]
    
    # do refine runs number of times
    for i in 1:runs
        # store all rotations against xyz axes
        rot_x = results_docking[1][1,:].R[1]
        rot_y = results_docking[1][1,:].R[2]
        rot_z = results_docking[1][1,:].R[3]

        # get random new numbers for rotations
        α = deg2rad(rand(Uniform(rot_x-Float32(20),rot_x+Float32(20))))
        β = deg2rad(rand(Uniform(rot_y-Float32(20),rot_y+Float32(20))))
        γ = deg2rad(rand(Uniform(rot_z-Float32(20),rot_z+Float32(20))))

        # create rotation 
        R = RotXYZ{Float32}(α,β,γ)
        
        # generate record and use last added rotation 
        # maybe nice
        record = generate_record(grid_A, R, roomcoordiantes_B, centroids, gridsize)
        
        # check if record is better than first one in
        # current results
        if (results_docking[1][1,:].score[1] < record.score)
            push!(results_docking[1], record)
            # sort scoring_table
            sort!(results_docking[1], [:score], rev=[true])
        end
    end
    
    # return new scoring table
    return results_docking[1]
end