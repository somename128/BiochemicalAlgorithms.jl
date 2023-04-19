using BiochemicalAlgorithms
using Meshes
using BenchmarkTools

include("create_atomballs.jl")
include("create_centroids.jl")

function grid_representation(protein)
    # calculate atomballs around proteins atoms
    atomballs = create_atomballs(protein)
    println(length(atomballs))
    # calculate centroids in a 128x128x128 grid with cells
    # of 1 angström 
    centroids = create_centroids(128,1)

    #extract min max (in rounded int) of atom coordinates of protein
    min_max = min_max_atoms(protein)
    min_x = min_max[1]
    max_x = min_max[2]
    min_y = min_max[3]
    max_y = min_max[4]
    min_z = min_max[5]
    max_z = min_max[6]

    println("Set marked cells...")
    # inside-outside
    # TODO: for loop with new colored_cells structure
    colored_cells = Base.Vector{Int64}()

    # check if centroids of cells inside balls of atoms
    # and store position of colored cell 
    @time begin
    for i in centroids, j in atomballs
        if(Base.in(i,j))
            # returns vector thats why position[1]
            # dont know if vector of vectors or number better 
            # for future calculations
            #
            # findall returns indice of i in centroids if a centroid i lies
            # in an atomball j -> stored in colored_cells if not already in storage
            position = findall(item -> item == i, centroids)
            if(!Base.in(position[1], colored_cells))
                push!(colored_cells,position[1])
            end
        end  
    end
    end

    println("Build 1D grid representation...")
    inner_outer_grid = zeros(N*N*N)

    # set the grid-position where "centroid inside atomball"
    # and all six cells around i are also into colored_cells (sourrounded by atoms)
    # to -15 (inside atom)
    # all other colored_cells are surface cells and set to 1
    # TODO: all surrounding cells
    for i in colored_cells
        if(Base.in(i-N*N, colored_cells) && Base.in(i-N,colored_cells) && Base.in(i-1, colored_cells)
            && Base.in(i+1, colored_cells) && Base.in(i+N, colored_cells) && Base.in(i+N*N, colored_cells))
            # inside
            inner_outer_grid[i] = -15
        else
            # surface
            inner_outer_grid[i] = 1
        end
    end

    # change 1D to 3D for FFTW library
    println("Build 3D grid representation...")
    inner_outer_grid_3D = zeros(N,N,N)

    for i in 1:N*N*N
        inner_outer_grid_3D[i] = inner_outer_grid[i]
    end

    return inner_outer_grid_3D
end