using BiochemicalAlgorithms
using Meshes

function create_atomballs(protein::Molecule{Float32})
    println("Extract room coordinates...")
    # extract room coordinates of atoms of the protein
    atoms_in_space = atoms_df(protein).r

    # transfer atom coordinates in mesh points
    atoms_in_space_points = Base.Vector{Meshes.Point3}()

    for i in atoms_in_space
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end

    # balls with radius r and atom points as center
    # TODO: different r for different atoms
    println("Build atomballs...")
    atomballs = Base.Vector{Meshes.Ball}() 
    r = 1.8

    for i in atoms_in_space_points
        b = Meshes.Ball(i, r)
        push!(atomballs, b)
    end

    return atomballs
end