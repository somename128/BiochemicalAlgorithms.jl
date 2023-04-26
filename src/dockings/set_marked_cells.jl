using BiochemicalAlgorithms

function set_marked_cells(atomballs,centroids,protein)
    # initialize vector with datatype of centroids
    colored_cells = Vector{Meshes.CartesianIndex{3}}()

    #extract min max (in rounded int) of atom coordinates of protein
    min_max = min_max_atoms(protein)
    min_x = min_max[1]
    max_x = min_max[2]
    min_y = min_max[3]
    max_y = min_max[4]
    min_z = min_max[5]
    max_z = min_max[6]

    # store centroids that are inside a atom radius in colored_cells
    for j in atomballs
        Threads.@threads for i in centroids[min_x:max_x,min_y:max_y,min_z:max_z]
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

    # return cells that represent protein in grid structure
    return colored_cells
end