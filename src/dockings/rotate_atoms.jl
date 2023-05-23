function rotate_atoms(roomcoordinates::Vector{Vector3{Float32}}, rotation::Matrix3{Float32})
    # needed because of race conditions of every thread
    # every thread makes a copy of roomcoordinates
    # otherwise error (but confused because reading ...)
    # maybe shared array
    atoms = roomcoordinates
    # function for broadcasting multiplication with rotation
    # matrix over every roomcoordinate
    f(v::Vector3{Float32}) = rotation * v 

    # use function
    atoms_rotated = f.(atoms)
    
    return atoms_rotated
end