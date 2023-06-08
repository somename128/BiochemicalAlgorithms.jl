using BiochemicalAlgorithms
using Meshes

include("extract_roomcoordinates.jl")

function create_atomballs(roomcoordinates::Vector{Vector3{Float32}})
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
    # for easy geometry for better couting of colored_cells 
    r = 0.6

    for i in atoms_in_space_points
        b = Meshes.Ball(i, r)
        push!(atomballs, b)
    end

    return atomballs
end