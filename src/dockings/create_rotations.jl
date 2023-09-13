using Quaternions

include("quaternion_functions.jl")

function create_rotations()
    
    # initialize arrays/vectors for rotation matrices and angles for rotations
    rotations = Vector{QuaternionF32}()
    angles_x = Array{Float32}(undef, 18)
    angles_y = Array{Float32}(undef, 10)
    angles_z = Array{Float32}(undef, 18)

    # fill array with angle values from 0 to 340 stepsize 20 (0==360)
    for i in 1:18
        a = deg2rad((i-1) * 20)
        angles_x[i] = a
        angles_z[i] = a
    end

    # fill array with angle values from 0 to 180 stepsize 20
    for i in 1:10
        angles_y[i] = deg2rad((i-1) * 20)
    end


    # loop over angles to store all possible rotations
    for γ in angles_z, β in angles_y, α in angles_x
        # α degrees rotation aroud x-axis
        q1 = quat_from_axisangle([1,0,0], α)
        # β degrees rotation aroud y-axis
        q2 = quat_from_axisangle([0,1,0], β)
        # β degrees rotation aroud z-axis
        q3 = quat_from_axisangle([0,0,1], γ)
        # combine quaternions for full rotation around xyz
        q = q1*q2*q3
        # add to array/vector of quaternions
        push!(rotations, q)
    end

    return rotations
end
