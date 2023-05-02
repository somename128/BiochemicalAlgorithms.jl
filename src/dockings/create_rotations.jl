using BiochemicalAlgorithms

function create_rotations()
    
    # initialize arrays/vectors for rigidtransforms and angles for transforms
    rigidtransforms = Vector{RigidTransform}()
    angles = Vector{Float32}()

    # fill array with angle values from 0 to 360 stepsize 5
    # for i in 1:73
    for i in 1:8
        a = (i-1) * 5
        push!(angles,a)
    end
    
    # translationvector (no translation needed)
    t = Vector3{Float32}(0,0,0)

    # loop over angles to store all possible rotations
    for α in angles, β in angles, γ in angles
        # rotation matrices per axe, R as resulting rotation around x y z
        R_x = Matrix3{Float32}([1 0 0; 0 cosd(α) -sind(α); 0 sind(α) cosd(α)])
        R_y = Matrix3{Float32}([cosd(β) 0 sind(β); 0 1 0; -sind(β) 0 cosd(β)])
        R_z = Matrix3{Float32}([cosd(γ) -sind(γ) 0; sind(γ) cosd(γ) 0; 0 0 1])
        R = R_x * R_y * R_z

        # create Biochemicals rigidtransform
        rigidtransform = RigidTransform(R,t)

        # add to array/vector of rigidtransformations
        push!(rigidtransforms, rigidtransform)
    end

    return rigidtransforms
end
