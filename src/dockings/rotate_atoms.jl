include("quaternion_functions.jl")

function rotate_atoms(roomcoordinates::Vector{Tuple{String, Vector3{Float32}}}, rotation::QuaternionF32, gridsize::Int32)
    # needed because of race conditions of every thread
    # every thread makes a copy of roomcoordinates
    # otherwise error (dont know why)
    # maybe shared array
    roomcoor = roomcoordinates
    # extract roomcoordinates from tuple
    atoms = Array{Vector3{Float32}}(undef, length(roomcoor))
    [atoms[i] = roomcoor[i][2] for i in eachindex(roomcoor)] 
    # idea:
    # 1. translate in origin (0,0,0) 
    # 2. rotate
    # 3. translate back in center of grid

    # 1. 
    # get mass center and broadcast translation to origin (0,0,0)
    mc = mass_center(atoms)
    for i in eachindex(atoms)
        atoms[i] = rotate_vector(rotation, (-1)*mc + atoms[i])
    end
    # f(v::Vector3{Float32}) = (-1)*mc + v
    # atoms = f.(atoms)
    # 2.
    # function for broadcasting rotation with quaternion over every 
    # roomcoordinate in origin
    # g(v::Vector3{Float32}) = rotate_vector(rotation, v)
    # atoms = g.(atoms)

    # 3.
    # get mass center and translate atoms back to center of grid
    mc = mass_center(atoms)
    #h(v::Vector3{Float32}) = (-1)*Vector3{Float32}(mc[1]-gridsize/2, mc[2]-gridsize/2, 
    #    mc[3]-gridsize/2) + v
    #atoms = h.(atoms)
    for i in eachindex(atoms)
        atoms[i] = (-1)*Vector3{Float32}(mc[1]-gridsize/2, mc[2]-gridsize/2, 
            mc[3]-gridsize/2) + atoms[i]
    end

    # put back together tuple
    rotated_tuple = [(roomcoor[i][1], atoms[i]) for i in eachindex(roomcoor)]
    
    return rotated_tuple
end