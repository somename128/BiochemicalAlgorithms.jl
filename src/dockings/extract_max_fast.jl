# function that extracts the maximum value of a complex or real 3D matrix
function extract_max_fast(C::Array{ComplexF32, 3})
    
    # return max value and the shifts 
    return Int32(argmax(real(C))[1]), Int32(argmax(real(C))[2]), Int32(argmax(real(C))[3]), Float32(maximum(real(C)))
end