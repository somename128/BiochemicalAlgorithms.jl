using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using ProgressBars
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations

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

N = Int32(64)
rotations = create_rotations()
R = Vector{Matrix3{Float32}}()
r = RotXYZ(deg2rad(-100),deg2rad(80),deg2rad(0))
push!(R,r)
centroids = create_centroids(N, one(Int32))
protein_A = load_and_trans_pdb("dummy_protein_vol2.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
protein_B = load_and_trans_pdb("dummy_ligand_vol2.pdb", N)
roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
A = grid_representation(roomcoordiantes_atoms_A, N, centroids)
B = grid_representation(roomcoordiantes_atoms_B, N, centroids)
shift = CartesianIndex(-61+N/2,-61+N/2,-3)
B_r = rotate_atoms(roomcoordiantes_atoms_B, R[1], N)
B_grid = grid_representation(B_r, N, centroids)
mcA = mass_center(roomcoordiantes_atoms_A)
mcB = mass_center(B_r)
# g(v::Vector3{Float32}) = (-1)*Vector3{Float32}(mcB[1]-33, mcB[2]-30, mcB[3]-2) + v
h(v::Vector3{Float32}) = (-1)*Vector3{Float32}(mcB[1]-N/2+2, mcB[2]-N/2+2, mcB[3]-N/2+2) + v
# atoms_translated_A = g.(roomcoordiantes_atoms_A)
atoms_translated_B = h.(B_r)

atoms_in_space_points = Base.Vector{Meshes.Point3}()
#=
# show grid rep
for i in CartesianIndices(A)
    if (A[i] != 0 && A[i] != -15)
        v = Meshes.Point(i[1],i[2],i[3])
        push!(atoms_in_space_points, v)
    end
end
=#
for i in CartesianIndices(B_grid)
    if (B_grid[i] != 0)
        v = Meshes.Point(i[1],i[2],i[3])
        #if (Base.in(v, atoms_in_space_points))
            push!(atoms_in_space_points, v)
        #end
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
    if (Base.in(v, atoms_in_space_points))
        push!(atoms_in_space_points, v)
    end
end
=#

viz(atoms_in_space_points, color = 1:length(atoms_in_space_points))
# viz(atoms_in_space_points)


