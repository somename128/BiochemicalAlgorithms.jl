using Meshes

function create_atomballs(roomcoordinates::Vector{Vector3{Float32}}, radius::Float32)
    # println("Extract room coordinates...")

    # transfer atom coordinates in mesh points
    atoms_in_space_points = Base.Vector{Meshes.Point3}()

    for i in roomcoordinates
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end

    # balls with radius r and atom points as center
    # TODO: different r for different atoms
    # println("Build atomballs...")
    atomballs = Base.Vector{Meshes.Ball}()

    # for easy geometry radius = 0.7 for better couting of colored_cells 
    # normally 1.8 Å or indiviual for every element
    for i in atoms_in_space_points
        b = Meshes.Ball(i, radius)
        push!(atomballs, b)
    end

    return atomballs
end