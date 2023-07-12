include("quaternion_functions.jl")

function rotate_atoms(roomcoordinates::Vector{Vector3{Float32}}, rotation::QuaternionF32, gridsize::Int32)
    # needed because of race conditions of every thread
    # every thread makes a copy of roomcoordinates
    # otherwise error (but confused because reading ...)
    # maybe shared array
    atoms = roomcoordinates
    # rotate atoms around xyz via quaternion
    g(v::Vector3{Float32}) = rotate_vector(rotation, v)
    atoms_rotated = g.(atoms)
    
    return atoms_rotated
end