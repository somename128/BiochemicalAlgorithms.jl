function create_inner_outer_grid(colored_cells::Vector{Int64}, N::Int64)  
    # println("Build 1D grid representation...")
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
    # println("Build 3D grid representation...")
    inner_outer_grid_3D = zeros(N,N,N)

    for i in 1:N*N*N
        inner_outer_grid_3D[i] = inner_outer_grid[i]
    end

    return inner_outer_grid_3D
end