using Combinatorics
using Quaternions
using Permutations
using DataFrames

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

    return (x-1)/res
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
    # encode all entries of v in a dictionary
    d = Dict{Int32, Float32}(1 => v[1], 2 => v[2], 3 => v[3], 4 => v[4])

    # new "encoded" vector for permutations
    w = Int32[1, 2, 3, 4]
    
    for i in permutations(w)
        # only push if its an even permutation
        #println(i) 
        #println(sign(Permutation(i)))
        if sign(Permutation(i)) == 1
            push!(array, QuaternionF32(d[i[1]], d[i[2]], d[i[3]], d[i[4]]))
        end 
    end
end


