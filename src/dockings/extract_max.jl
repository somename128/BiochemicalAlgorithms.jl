using TypedTables
using BenchmarkTools

include("helpers.jl")

# function that extracts the maximum value of a complex or real 3D matrix
function extract_max(C)
    C_real = real(C)
    
    # initialize record for max
    # to get access to first element filled with one element
    t = Table(α=[0.0], β= [0.0], γ=[0.0], score=[0.0])

    C_sparse = zero_small!(C_real,0.5)

    # walk through all values in 3D matrix and store them in table
    # indices and score
    for i in CartesianIndices(C_sparse)
        if (t[1].score < C_sparse[i])
            t[1] = (α=i[1], β=i[2], γ=i[3], score=C_sparse[i])
        end
    end

    # find and store maximum of all scores
    # for now only first max element stored
    # max = t[findmax(t.score)[2]]

    # return max value and the shifts 
    return t[1]
end