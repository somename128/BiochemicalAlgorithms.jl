using Combinatorics
using Quaternions

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
# if value is greater than gridsize*res/2
# values are interpreted as negativ shifts (value - gridsize*res)
# to transform in 1 angstöm steps for using in angstöm space
# x is returned as division through resolution value
function interp!(x::Int32, N::Int32, res::Int32)
    if x > N*res/2
        x -= N*res
    end

    return x/res
end

# function that takes a vector, calculates all permutations
# convert the result to quaternions and push them in the given array
function perm_quat!(v::Vector{Float32}, array::Vector{QuaternionF32})
    for i in permutations(v)
        push!(array, QuaternionF32(i[1], i[2], i[3], i[4]))
    end
end

# function like above, but only pushes even permutations
function perm_quat_even!(v::Vector{Float32}, array::Vector{QuaternionF32})
    for i in permutations(v)
        errors = 0
        for j in 1:4
            if i[j] != v[j]
                errors += 1
            end
        end

        # only push if its an even permutation (free interpretation
        # of even)
        if (-1)^errors == 1
            push!(array, QuaternionF32(i[1], i[2], i[3], i[4]))
        end 
    end
end