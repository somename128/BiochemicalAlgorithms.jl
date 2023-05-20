using BiochemicalAlgorithms

function create_rotations()
    
    # initialize arrays/vectors for rigidtransforms and angles for transforms
    rigidtransforms = Vector{RigidTransform{Float32}}()
    angles_x = Vector{Float32}()
    angles_y = Vector{Float32}()
    angles_z = Vector{Float32}()

    # fill array with angle values from 0 to 360 stepsize 5
    #for i in 1:73
    for i in 1:10
        a = (i-1) * 5
        push!(angles_x,a)
        push!(angles_z,a)
    end

    # fill array with angle values from 0 to 180 stepsize 5
    # for i in 1:37
    for i in 1:10
        a = (i-1) * 5
        push!(angles_y,a)
    end

    
    # translationvector (no translation needed)
    t = Vector3{Float32}(0,0,0)

    # loop over angles to store all possible rotations
    for ϕ in angles_z, Θ in angles_y, ψ in angles_x
        # rotation matrix after xyz convention (see Goldsteins Classical Mechanics(1980))
        R = Matrix3{Float32}([cosd(Θ)*cosd(ϕ) cosd(Θ)*sind(ϕ) -sind(Θ); 
            sind(ψ)*sind(Θ)*cosd(ϕ)-cosd(ψ)*sind(ϕ) sind(ψ)*sind(Θ)*sind(ϕ)+cosd(ψ)*cosd(ϕ) cosd(Θ)*sind(ψ); 
            cosd(ψ)*sind(Θ)*cos(ϕ)+sind(ψ)*sind(Θ) cosd(ψ)*sind(Θ)*sind(ϕ)-sind(ψ)*cosd(ϕ) cosd(Θ)*cosd(ψ)])

        # create Biochemicals rigidtransform
        rigidtransform = RigidTransform{Float32}(R,t)

        # add to array/vector of rigidtransformations
        push!(rigidtransforms, rigidtransform)
    end

    return rigidtransforms
end
