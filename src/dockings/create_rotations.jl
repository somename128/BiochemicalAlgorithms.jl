using BiochemicalAlgorithms
using Rotations

function create_rotations()
    
    # initialize arrays/vectors for rotation matrices and angles for rotations
    rotations = Vector{Matrix3{Float32}}()
    angles_x = Vector{Float32}()
    angles_y = Vector{Float32}()
    angles_z = Vector{Float32}()

    # fill array with angle values from 0 to 360 stepsize 20
    for i in 1:19
        a = (i-1) * 5
        push!(angles_x,a)
        push!(angles_z,a)
    end

    # fill array with angle values from 0 to 180 stepsize 20
    for i in 1:10
        a = (i-1) * 5
        push!(angles_y,a)
    end


    # loop over angles to store all possible rotations
    for ϕ in angles_z, Θ in angles_y, ψ in angles_x
        # rotation matrix after xyz convention (see Goldsteins Classical Mechanics(1980))
        R = RotYXZ(Θ,ψ,ϕ)
        # add to array/vector of rigidtransformations
        push!(rotations, R)
    end

    return rotations
end
