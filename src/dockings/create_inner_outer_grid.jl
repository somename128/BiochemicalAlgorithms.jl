function create_inner_outer_grid(colored_cells::Vector{CartesianIndex}, N::Int64)  
    # intialize 3D grid witz zeros
    inner_outer_grid = zeros(N,N,N)

    # set the grid-position where "centroid inside atomball"
    # and all six cells around i are also into colored_cells (sourrounded by atoms)
    # to -15 (inside atom)
    # all other colored_cells are surface cells and set to 1
    # TODO: all surrounding cells
    for i in colored_cells
        if(Base.in(CartesianIndex(i[1]-1,i[2],i[3]), colored_cells) 
            && Base.in(CartesianIndex(i[1]+1,i[2],i[3]),colored_cells) 
            && Base.in(CartesianIndex(i[1],i[2]-1,i[3]), colored_cells)
            && Base.in(CartesianIndex(i[1],i[2]+1,i[3]), colored_cells) 
            && Base.in(CartesianIndex(i[1],i[2],i[3]-1), colored_cells) 
            && Base.in(CartesianIndex(i[1],i[2],i[3]+1), colored_cells))
            # inside
            inner_outer_grid[i] = -15
        else
            # surface
            inner_outer_grid[i] = 1
        end
    end

    return inner_outer_grid
end