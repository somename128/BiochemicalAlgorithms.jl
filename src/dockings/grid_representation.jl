using BiochemicalAlgorithms
using Meshes
using BenchmarkTools

function grid_representation(protein)
    println("Extract room coordinates...")
    # extract room coordinates of atoms of the protein
    atoms_in_space = protein.atoms.r

    # transfer atom coordinates in mesh points
    atoms_in_space_points = Base.Vector{Meshes.Point3}()

    for i in atoms_in_space
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end

    println("Build grid...")
    # grid properties
    lower_left = (0,0,0)
    N = 128
    spcing_factor = 1
    upper_right = (N,N,N) 
    spcing = (N*spcing_factor, N*spcing_factor, N*spcing_factor)

    # create 3D grid with N/N*spcing_factor spacing in each dimension, origin at (0,0,0)
    grid = Meshes.CartesianGrid(lower_left, upper_right, dims=spcing)
    # centroid of each cell
    println("Build centroids...")
    centroids = Meshes.centroid.(grid)

    # balls with radius r and atom points as center
    # TODO: different r for different atoms
    println("Build atomballs...")
    atomballs = Base.Vector{Meshes.Ball}()
    r = 1.8

    for i in atoms_in_space_points
        b = Meshes.Ball(i, r)
        push!(atomballs, b)
    end

    @time begin
    println("Set marked cells...")
    #inside-outside
    colored_cells = Base.Vector{Int64}()

    # check if centroids of cells inside balls of atoms
    # and store position of colored cell 
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