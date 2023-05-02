# function to set all very small values to zero
function zero_small!(M, tol)
    for ι in eachindex(M)
        if abs(M[ι]) ≤ tol
            M[ι] = 0
        end
    end
    M
end

# functions that check if protein boundaries are out of bounds
function less_than_one!(x::Int64)
    if x < 1
        x = 1
    end

    return x
end
