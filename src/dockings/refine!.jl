using Rotations
using Distributions

include("generate_record.jl")
include("quaternion_functions.jl")

function refine!(results_docking::Tuple, runs::Int32)

    # store grid of protein A, roomcoordinates of protein b,
    # centroids and gridsize
    grid_A = results_docking[2]
    roomcoordiantes_B = results_docking[3]
    centroids = results_docking[4]
    gridsize = results_docking[5]
    
    # do refine runs number of times
    for i in 1:runs
        # store all rotations against xyz axes in quaternions
        # and generate quaternion for rotation
        q_x = quat_from_axisangle([1,0,0], deg2rad(results_docking[1][1,:].R[1]))
        q_y = quat_from_axisangle([0,1,0], deg2rad(results_docking[1][1,:].R[2]))
        q_z = quat_from_axisangle([0,0,1], deg2rad(results_docking[1][1,:].R[3]))
        q = q_x*q_y*q_z

        # get axis and angle from the best rotation (represented in quaternion)
        quat_axis = axisangle_from_quat(q)[1]
        quat_angle = axisangle_from_quat(q)[2]

        # get new quat_axis depending on old one?
        # Kent or Bingham dirstribution for new quat_axis?
        # new_quat_axis
        # new_quat_angle

        # create rotation 
        R = quat_from_axisangle(quat_axis, quat_angle)
        
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