using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools
using DelimitedFiles


@time begin
# function for calculating grid representation
function grid_representation(path_to_pdb, translation_vector::Vector3{Float32})
    # load protein data from PDB
    protein = load_pdb(path_to_pdb)
    println("PDB loaded...")

    # translation vector
    t = translation_vector

    # translate protein in "positive space"
    # TODO: better choice of translation vector
    BiochemicalAlgorithms.translate!(protein,t)

    # extract room coordinates of atoms of the protein
    atoms_in_space = protein.atoms.r

    # transfer atom coordinates in mesh points
    atoms_in_space_points = Base.Vector{Meshes.Point3}()

    for i in atoms_in_space
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end

    # visualize atoms of protein
    # viz(atoms_in_space_points)

    #=
    # for the start: hard coded grid and translated protein
    # TODO: improvement

    # initalize vectors for storing x y z coordinates seperately
    X = Vector{Float32}()
    Y = Vector{Float32}()
    Z = Vector{Float32}()

    # fill vectors with coordinates
    for i in atoms_in_space
        push!(X, i[1])
        push!(Y, i[2])
        push!(Z, i[3])
    end

    # calculate x y z min and max for grid size
    min_x = minimum(X)
    max_x = maximum(X)
    min_y = minimum(Y)
    max_y = maximum(Y)
    min_z = minimum(Z)
    max_z = maximum(Z)


    lower_left = Meshes.Point(max_x, min_y, min_z)
    upper_right = Meshes.Point(min_x, max_y, max_z)
    spacing = (100,100,100)
    =#

    # grid properties
    lower_left = (0,0,0)
    N = 64
    spcing_factor = 1
    upper_right = (N,N,N) 
    spcing = (N*spcing_factor, N*spcing_factor, N*spcing_factor)

    # create 3D grid with N/N*spcing_factor spacing in each dimension, origin at (0,0,0)
    grid = Meshes.CartesianGrid(lower_left, upper_right, dims=spcing)
    println(string(N) * "x" * string(N) * "x" * string(N) * " Grid built...")
    # centroid of each cell
    centroids = Meshes.centroid.(grid)
    println("Centroids loaded...")

    # balls with radius r and atom points as center
    # TODO: different r for different atoms 
    atomballs = Base.Vector{Meshes.Ball}()
    r = 1.8

    for i in atoms_in_space_points
        b = Meshes.Ball(i, r)
        push!(atomballs, b)
    end
    println(string(length(atomballs)) * " Atomballs built...")

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
    println(string(length(colored_cells)) * " Colored cells found...")

    inner_outer_grid = zeros(N*N*N)

    # set the grid-position where "centroid inside atomball"
    # and all six cells around i are also into colored_cells (sourrounded by atoms)
    # to -15 (inside atom)
    # all other colored_cells are surface cells and set to 1
    for i in colored_cells
        if(Base.in(i-N*N, colored_cells) && Base.in(i-N,colored_cells) && Base.in(i-1, colored_cells)
            && Base.in(i+1, colored_cells) && Base.in(i+N, colored_cells) && Base.in(i+N*N, colored_cells))
            inner_outer_grid[i] = -15
            # println("inside!")
        else
            inner_outer_grid[i] = 1
            # println("surface!")
        end
    end
    println("Inner-outer-Grid built...")

    return inner_outer_grid
end

t1 = Vector3{Float32}(20,0,0)
t2 = Vector3{Float32}(0,-40,0)

A = grid_representation("2ptc_protein.pdb", t1)
B = grid_representation("2ptc_ligand.pdb", t2)

println("Start writing grid A and B in txt files...")
# write gird representation in txt file
writedlm("grid_A.txt", A)
writedlm("grid_B.txt", B)

println("Start correlation...")
# correlation function
C = dot(A,B)

println("Scoring: " * string(C))

end