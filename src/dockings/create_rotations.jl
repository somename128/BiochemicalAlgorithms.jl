using Quaternions

include("quaternion_functions.jl")

function create_rotations()
    
    # initialize arrays/vectors for rotation matrices and angles for rotations
    rotations = Vector{QuaternionF32}()
    angles_x = Vector{Float32}()
    angles_y = Vector{Float32}()
    angles_z = Vector{Float32}()

    # fill array with angle values from 0 to 360 stepsize 20
    for i in 1:19
        a = deg2rad((i-1) * 20)
        push!(angles_x,a)
        push!(angles_z,a)
    end

    # fill array with angle values from 0 to 180 stepsize 20
    for i in 1:10
        a = deg2rad((i-1) * 20)
        push!(angles_y,a)
    end


    # loop over angles to store all possible rotations
    for γ in angles_z, β in angles_y, α in angles_x
        # α degrees rotation aroud x-axis
        q1 = quat_from_axisangle([1,0,0], deg2rad(α))
        # β degrees rotation aroud y-axis
        q2 = quat_from_axisangle([0,1,0], deg2rad(β))
        # β degrees rotation aroud z-axis
        q3 = quat_from_axisangle([0,0,1], deg2rad(γ))
        # combine quaternions for full rotation around xyz
        q = q1*q2*q3
        # add to array/vector of quaternions
        push!(rotations, q)
    end

    return rotations
end
