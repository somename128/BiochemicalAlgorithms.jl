using DataFrames
using BenchmarkTools

include("helpers.jl")

# function that extracts the maximum value of a complex or real 3D matrix
function extract_max(C)
    C_real = real(C)
    
    # initialize record for max
    # to get access to first element filled with one element
    t = DataFrame(α=zero(Int32), β=zero(Int32), γ=zero(Int32), score=zero(Float32))

    C_sparse = zero_small!(C_real,Float32(0.5))

    # walk through all values in 3D matrix and store them in table
    # indices and score
    for i in CartesianIndices(C_sparse)
        if (t.score[1] < C_sparse[i])
            t[1,:] = (α=Int32(i[1]), β=Int32(i[2]), γ=Int32(i[3]), score=Float32(C_sparse[i]))
        end
    end

    # find and store maximum of all scores
    # for now only first max element stored
    # max = t[findmax(t.score)[2]]

    # return max value and the shifts 
    return t
end