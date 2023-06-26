using FFTW
using BiochemicalAlgorithms
using DSP
using FourierTools
using DataFrames
using SparseArrays

include("create_rotations.jl")
include("get_degrees.jl")
include("generate_record.jl")

N = 4

A = zeros(Float32, N,N,N)
B = zeros(Float32, N,N,N)

for i in CartesianIndices(A[2:3,2:3,2:3])
    A[i+CartesianIndex(1,1,1)] = 1
end

for i in eachindex(B)
    if (i%2==0)
        B[i] = -i
    else
        B[i] = i
    end
end

C = ifft(fft(A).*fft(B))

max = extract_max(C)
println(max.α[1]-N," ",max.β[1]-N," ",max.γ[1]-N)

