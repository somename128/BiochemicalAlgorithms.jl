using BiochemicalAlgorithms

function extract_roomcoordinates(protein::Molecule{Float32})
    # store roomcoordinates
    atoms_r = atoms_df(protein).r
    # initialize vector for roomcoordinates of atoms
    roomcoordinates_atoms = Vector{Vector3{Float32}}()

    # extract the roomcoordinates from dataframe
    for i in eachindex(atoms_r)
        push!(roomcoordinates_atoms, atoms_r[i])
    end

    return roomcoordinates_atoms

end