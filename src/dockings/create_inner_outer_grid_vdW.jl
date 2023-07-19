function create_inner_outer_grid_vdW(inner_cells::Vector{Int32}, surface_cells::Vector{Int32}, N::Int32, res::Int32)
    # initalize grid in 1D
    # size needed is M = N*res
    M = N*res
    inner_outer_grid = zeros(Float32, M*M*M)

    # set inner grid cells to -15
    for i in inner_cells
        inner_outer_grid[i] = Float32(-15)
    end

    # set surface grid cells to 1
    for i in surface_cells
        inner_outer_grid[i] = one(Float32)
    end

    # change 1D to 3D for FFTW library
    inner_outer_grid_3D = zeros(Float32, M,M,M)

    for i in 1:M*M*M
        inner_outer_grid_3D[i] = inner_outer_grid[i]
    end

    return inner_outer_grid_3D
end

