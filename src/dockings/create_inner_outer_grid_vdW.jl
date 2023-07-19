function create_inner_outer_grid_vdW(inner_cells::Vector{Int32}, surface_cells::Vector{Int32}, N::Int32)
    # initalize grid in 1D
    inner_outer_grid = zeros(Float32, N*N*N)

    # set inner grid cells to -15
    for i in inner_cells
        inner_outer_grid[i] = Float32(-15)
    end

    # set surface grid cells to 1
    for i in surface_cells
        inner_outer_grid[i] = one(Float32)
    end

    # change 1D to 3D for FFTW library
    inner_outer_grid_3D = zeros(Float32, N,N,N)

    for i in 1:N*N*N
        inner_outer_grid_3D[i] = inner_outer_grid[i]
    end

    return inner_outer_grid_3D
end

