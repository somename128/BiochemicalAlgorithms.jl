function rotate_atoms(roomcoordinates::Vector{Vector3{Float32}}, rotation::Matrix3{Float32})
    atoms = roomcoordinates
    atoms_rotated = Vector{Vector3{Float32}}()

    for i in eachindex(atoms)
        v = rotation * atoms[i]
        push!(atoms_rotated,v)
    end
    
    return atoms_rotated
end