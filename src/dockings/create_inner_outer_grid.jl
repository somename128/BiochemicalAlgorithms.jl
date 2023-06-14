function create_inner_outer_grid(colored_cells::Vector{Int32}, N::Int32)  
    # println("Build 1D grid representation...")
    inner_outer_grid = zeros(Float32, N*N*N)

    # set the grid-position where "centroid inside atomball"
    # and all six cells around i are also into colored_cells (surrounded by atoms)
    # to -15 (inside atom)
    # all other colored_cells are surface cells and set to 1
    for i in colored_cells
        if (# under front left (-z,-y,-x direction)
            Base.in(i-N*N-N-1, colored_cells)
            # under front (-z,-y direction)
            && Base.in(i-N*N-N, colored_cells)
            # under front right (-z,-y,+x direction)
            && Base.in(i-N*N-N+1, colored_cells)
            # under left (-z,-x direction)
            && Base.in(i-N*N-1, colored_cells)
            # under (-z direction)
            && Base.in(i-N*N, colored_cells)
            # under right (-z,+x direction)
            && Base.in(i-N*N+1, colored_cells)
            # under behind left (-z,+y,-x direction)
            && Base.in(i-N*N+N-1, colored_cells)
            # under behind (-z,+y direction)
            && Base.in(i-N*N+N, colored_cells)
            # under behind right (-z,+y,+x direction)
            && Base.in(i-N*N+N+1, colored_cells)
            # front left (-y,-x direction)
            && Base.in(i-N-1,colored_cells)
            # front (-y direction)
            && Base.in(i-N,colored_cells)
            # front right (-y,+x direction)
            && Base.in(i-N+1,colored_cells)
            # left (-x direction)
            && Base.in(i-1, colored_cells)
            # right (+x direction)
            && Base.in(i+1, colored_cells)
            # behind left (+y,-x direction)
            && Base.in(i+N-1, colored_cells)
            # behind (+y direction)
            && Base.in(i+N, colored_cells)
            # behind right (+y,+x direction)
            && Base.in(i+N+1, colored_cells)
            # over front left (+z,-y,-x direction)
            && Base.in(i+N*N-N-1, colored_cells)
            # over front (+z,-y direction)
            && Base.in(i+N*N-N, colored_cells)
            # over front right (+z,-y,+x direction)
            && Base.in(i+N*N-N+1, colored_cells)
            # over left (+z,-x direction)
            && Base.in(i+N*N-1, colored_cells)
            # over (+z direction)
            && Base.in(i+N*N, colored_cells)
            # over right (+z,+x direction)
            && Base.in(i+N*N+1, colored_cells)
            # over behind left (+z,+y,-x direction)
            && Base.in(i+N*N+N-1, colored_cells)
            # over behind (+z,+y direction)
            && Base.in(i+N*N+N, colored_cells)
            # over behind right (+z,+y,+x direction)
            && Base.in(i+N*N+N+1, colored_cells))
            # inside
            inner_outer_grid[i] = Float32(-15)
        else
            # surface
            inner_outer_grid[i] = one(Float32)
        end
    end

    # change 1D to 3D for FFTW library
    # println("Build 3D grid representation...")
    inner_outer_grid_3D = zeros(Float32, N,N,N)

    for i in 1:N*N*N
        inner_outer_grid_3D[i] = inner_outer_grid[i]
    end

    return inner_outer_grid_3D
end