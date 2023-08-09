using BiochemicalAlgorithms

function extract_roomcoordinates(protein::Molecule{Float32})
    # store roomcoordinates
    atoms = atoms_df(protein)

    # extract elements and roomcoordinates from dataframe
    roomcoordinates_atoms = [(i.name, i.r) for i in eachrow(atoms)]

    return roomcoordinates_atoms

end