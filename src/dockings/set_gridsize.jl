using BiochemicalAlgorithms

include("extract_roomcoordinates.jl")
include("extract_min_max.jl")

function set_gridsize(path_protein_A::String, path_protein_B::String)
    # extract roomcoordinates of atoms
    protein_A = molecules(load_pdb(path_protein_A))[1]
    protein_B = molecules(load_pdb(path_protein_B))[1]
    atoms_A = extract_roomcoordinates(protein_A)
    atoms_B = extract_roomcoordinates(protein_B)
    # extract roomcoordinates from tuples
    atomsA = Vector{Vector3{Float32}}()
    [push!(atomsA, i[2]) for i in atoms_A]
    atomsB = Vector{Vector3{Float32}}()
    [push!(atomsB, i[2]) for i in atoms_B]

    # function for extracting min and max in xyz direction
    min_max_A = extract_min_max(atomsA)
    min_max_B = extract_min_max(atomsB)

    # get the diameters via maximum - minimum for each dimension xyz
    diameters = Vector{Float32}()
    push!(diameters, min_max_A[2] - min_max_A[1])
    push!(diameters, min_max_A[4] - min_max_A[3])
    push!(diameters, min_max_A[6] - min_max_A[5])
    push!(diameters, min_max_B[2] - min_max_B[1])
    push!(diameters, min_max_B[4] - min_max_B[3])
    push!(diameters, min_max_B[6] - min_max_B[5])
    
    # store the greatest diameter plus buffer
    max_d = maximum(diameters) + 5

    # set the gridsize depending on the double value of the greatest
    # diameter (for fft)
    gridsize = zero(Int32)

    if (2*max_d <= 128 && 2*max_d > 64)
        gridsize = Int32(128)
    elseif (2*max_d <= 64 && 2*max_d > 32)
        gridsize = Int32(64)
    elseif (2*max_d <= 32 && 2*max_d > 16)
        gridsize = Int32(32)
    elseif (2*max_d <= 16 && 2*max_d > 8)
        gridsize = Int32(16)
    elseif (2*max_d <= 8 && 2*max_d > 4)
        gridsize = Int32(8)
    elseif (2*max_d < 4)
        gridsize = Int32(4)    
    else
        println("Your proteins are to heavy for my program.")
    end
    
    return gridsize
end