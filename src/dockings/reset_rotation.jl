using BiochemicalAlgorithms

function reset_rotation!(protein, rotation)
    # set translation vector for rigidtransform 
    # (not needed only for rotation)
    t = Vector3{Float32}(0,0,0)
    # create transposed (reverse) rotation matrix for resetting
    # rotation
    R = transpose(rotation.rotation)
    # set rigidtransform for reverse rotation
    rotation = RigidTransform{Float32}(R,t)
    # reset this steps protein rotation
    rigid_transform!(protein,rotation)
end