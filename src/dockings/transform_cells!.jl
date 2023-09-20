function transform_cells!(grid::Array{ComplexF32, 3})
    # cube that is shifted over the other cube
    kernel = ones(Float32, 3, 3, 3)
    tmp = Array{Float32, 3}(undef, 3, 3, 3)

    for i in CartesianIndices(grid[2:size(grid,1)-1, 2:size(grid,2)-1, 2:size(grid,3)-1])
        if (grid[i] == 1)
            # set small piece of greater grid
            for j in CartesianIndices(tmp)
                tmp[j] = real(grid[i+j-CartesianIndex(2,2,2)])
            end
            # calculate convolution            
            conv = sum(tmp .* kernel)
            # if all surrounding cells are one set cell to inner cell value
            # -15
            if (conv == 27)
                grid[i] = ComplexF32(-15)
            end
        end
    end

    return grid
end