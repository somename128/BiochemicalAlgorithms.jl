using FFTW
using BiochemicalAlgorithms
using DSP
using FourierTools
using DataFrames
using SparseArrays
using Meshes, MeshViz
using Makie, WGLMakie

include("create_rotations.jl")
include("get_degrees.jl")
include("generate_record.jl")
include("helpers.jl")

N = 2

A = zeros(Float32, N,N,N)
B = zeros(Float32, N,N,N)

for i in eachindex(A)
    A[i] = i
end
#=
for i in CartesianIndices(B[3:4,3:4])
    B[i+CartesianIndex(2,2)] = 8
end
=#
B[2,2,2] = 1
atoms_in_space_points = Base.Vector{Meshes.Point3}()
# show grid rep

for i in CartesianIndices(A)
    if (A[i] != 0)
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end
end
shift = CartesianIndex(-1,-1,-1)
for i in CartesianIndices(B)
    if (B[i] != 0)
        v = Meshes.Point(i[1]+shift[1],i[2]+shift[2],i[3]+shift[3])
        push!(atoms_in_space_points, v)
    end
end
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
# display(zero_small!(real(C),0.5))

viz(atoms_in_space_points, color = 1:length(atoms_in_space_points))
# viz(atoms_in_space_points)

display(A)
