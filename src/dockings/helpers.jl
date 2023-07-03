# function to set all very small values to zero
function zero_small!(M, tol::Float32)
    for ι in eachindex(M)
        if abs(M[ι]) ≤ tol
            M[ι] = 0
        end
    end
    M
end

# function to set a value to zero if its almost zero
function zero_small_one_value(x, tol::Float32)
    if abs(x) ≤ tol
        return 0
    end
    
    return x
end


# functions that check if protein boundaries are out of bounds
function less_than_one!(x::Int64)
    if x < 1
        x = 1
    end

    return x
end

# function that interprets value from fft
# if value is greater than gridsize/2
# values are interpret as negativ shifts 
function interp!(x::Int32, N::Int32)
    if x > N/2
        x -= N
    end

    return x
end