using FFTW
using BiochemicalAlgorithms
using DSP
using FourierTools
using DataFrames
using SparseArrays

include("create_rotations.jl")
include("get_degrees.jl")
include("generate_record.jl")
include("helpers.jl")

N = 2

A = zeros(Float32, N,N)
B = zeros(Float32, N,N)
#=
for i in CartesianIndices(A[2:3,1:2])
    A[i+CartesianIndex(1,0)] = 1
end

for i in CartesianIndices(B[3:4,3:4])
    B[i+CartesianIndex(2,2)] = 8
end
=#
A[1,1] = 8
B[1,2] = 1
#=
min_x = 0
min_y = 0
min_z = 0
# find first element
for i in CartesianIndices(A)
    if (A[i] != 0)
        global min_x = i[1]
        global min_y = i[2]
        global min_z = i[3]
        break
    end
end
=#


C = ifft(fft(A).*fft(B))

maxi = extract_max(C)
# println((1+maxi.α[1])%N," ",(2+maxi.β[1])%N," ",(1+maxi.γ[1])%N)
display(zero_small!(real(C),0.5))
