using SparseArrays

function create_inner_outer_grid(colored_cells::Vector{Int32}, N::Int32, is_smaller::Bool)  
    # println("Build 1D grid representation...")
    inner_outer_grid = zeros(Float32, N*N*N)

    # destingish between protein and ligand
    # if ligand set all colored cells to one 
    # the protein gets the destingtion between inner and surface
    if (!is_smaller)        
        # set the grid-position where "centroid inside atomball"
        # and all six cells around i are also into colored_cells (surrounded by atoms)
        # to -15 (inside atom) (for protein)
        # all other colored_cells are surface cells and set to 1
        for i in colored_cells
            if (# upper front left (-z,-y,-x direction)
                Base.in(i-N*N-N-1, colored_cells)
                # front left (-z,-y direction)
                && Base.in(i-N*N-N, colored_cells)
                # lower front left (-z,-y,+x direction)
                && Base.in(i-N*N-N+1, colored_cells)
                # upper front (-z,-x direction)
                && Base.in(i-N*N-1, colored_cells)
                # front (-z direction)
                && Base.in(i-N*N, colored_cells)
                # lower front (-z,+x direction)
                && Base.in(i-N*N+1, colored_cells)
                # upper front right (-z,+y,-x direction)
                && Base.in(i-N*N+N-1, colored_cells)
                # front right (-z,+y direction)
                && Base.in(i-N*N+N, colored_cells)
                # lower front right (-z,+y,+x direction)
                && Base.in(i-N*N+N+1, colored_cells)
                # upper left (-y,-x direction)
                && Base.in(i-N-1,colored_cells)
                # left (-y direction)
                && Base.in(i-N,colored_cells)
                # lower left (-y,+x direction)
                && Base.in(i-N+1,colored_cells)
                # upper (-x direction)
                && Base.in(i-1, colored_cells)
                # under (+x direction)
                && Base.in(i+1, colored_cells)
                # upper right (+y,-x direction)
                && Base.in(i+N-1, colored_cells)
                # right (+y direction)
                && Base.in(i+N, colored_cells)
                # lower right (+y,+x direction)
                && Base.in(i+N+1, colored_cells)
                # upper back left (+z,-y,-x direction)
                && Base.in(i+N*N-N-1, colored_cells)
                # back left (+z,-y direction)
                && Base.in(i+N*N-N, colored_cells)
                # lower back left (+z,-y,+x direction)
                && Base.in(i+N*N-N+1, colored_cells)
                # upper back (+z,-x direction)
                && Base.in(i+N*N-1, colored_cells)
                # back (+z direction)
                && Base.in(i+N*N, colored_cells)
                # lower back (+z,+x direction)
                && Base.in(i+N*N+1, colored_cells)
                # upper back right (+z,+y,-x direction)
                && Base.in(i+N*N+N-1, colored_cells)
                # back right (+z,+y direction)
                && Base.in(i+N*N+N, colored_cells)
                # lower back right (+z,+y,+x direction)
                && Base.in(i+N*N+N+1, colored_cells))
                # inside
                inner_outer_grid[i] = Float32(-15)
            else
                # surface
                inner_outer_grid[i] = one(Float32)
            end
        end
    else
        for i in colored_cells
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