using Rotations

function create_rotations()
    
    # initialize arrays/vectors for rotation matrices and angles for rotations
    rotations = Vector{RotXYZ{Float32}}()
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
        # rotation matrix after xyz convention (see Goldsteins Classical Mechanics(1980))
        R = RotXYZ{Float32}(α,β,γ)
        # add to array/vector of rigidtransformations
        push!(rotations, R)
    end

    return rotations
end
