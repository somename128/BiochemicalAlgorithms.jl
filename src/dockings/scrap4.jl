using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations
using DelimitedFiles

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("create_atomballs.jl")
include("correlation_docking.jl")
include("mass_center.jl")
include("rotate_atoms.jl")
include("helpers.jl")
include("extract_max.jl")

N = Int32(32)
rotations = create_rotations()
R = Vector{Matrix3{Float32}}()
r = RotXYZ(deg2rad(0),deg2rad(0),deg2rad(-140))
push!(R,r)
centroids = create_centroids(N, one(Int32))
protein_A = load_and_trans_pdb("dummy_protein.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
protein_B = load_and_trans_pdb("dummy_ligand.pdb", N)
roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
A = grid_representation(roomcoordiantes_atoms_A, N, centroids, false)
# B = grid_representation(roomcoordiantes_atoms_B, N, centroids, true)
# shift = CartesianIndex(-1, -1, -1)
B_r = rotate_atoms(roomcoordiantes_atoms_B, R[1], N)
h(v::Vector3{Float32}) = Vector3{Float32}(-1, -1, -7) + v
atoms_translated_B = h.(B_r)
B_grid = grid_representation(atoms_translated_B, N, centroids, true)

# scoring = Base.Vector{Meshes.Point3}()
atoms_in_space_points = Base.Vector{Meshes.Point3}()


# show grid rep
for i in CartesianIndices(A)
    if (A[i] != 0 && A[i] != -15)
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end
end

for i in CartesianIndices(B_grid)
    if (B_grid[i] != 0)
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end
end

#=
# way of rotations
for i in eachindex(rotations)
    B_r = rotate_atoms(roomcoordiantes_atoms_B,rotations[i], N)

    for a in B_r
        v = Meshes.Point(a[1],a[2],a[3])
        if (!Base.in(v, atoms_in_space_points))
            push!(atoms_in_space_points, v)
        end
    end
    println(i,"/",length(rotations))
end

=# 
#=
# show protein and rotated ligand
for i in roomcoordiantes_atoms_A
    v = Meshes.Point(i[1],i[2],i[3])
    push!(atoms_in_space_points, v)
end

for i in atoms_translated_B
    v = Meshes.Point(i[1],i[2],i[3])
    push!(atoms_in_space_points, v)
end
=#

# C = ifft(fft(A).*conj(fft(B_grid)))
#=
for i in CartesianIndices(C)
    v = Meshes.Point(i[1],i[2],real(C[i]))
    push!(scoring, v)
end
=#

# println(extract_max(zero_small!(real(C),Float32(0.5))))

viz(atoms_in_space_points, color = 1:length(atoms_in_space_points))
# viz(atoms_in_space_points)

#viz(scoring, color = 1:length(scoring))

