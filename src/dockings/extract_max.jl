using TypedTables

include("helpers.jl")

# function that extracts the maximum value of a complex 3D matrix
function extract_max(C)
    C_real = real(C)
    
    t = Table(α=[], β= [], γ=[], score=[])

    C_sparse = zero_small!(C_real,0.5)

    # walk through all values in 3D matrix and store them in table
    for k = 1:size(C_sparse,3)
        for j = 1:size(C_sparse,2)
            for i = 1:size(C_sparse,1)
                push!(t,(α=i, β=j, γ=k, score=C_sparse[i,j,k]))
            end
        end
    end

    # find and store maximum of all scores
    max = t[findmax(t.score)[2]]

   return max
end