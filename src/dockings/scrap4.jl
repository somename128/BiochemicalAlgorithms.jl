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
centroids = create_centroids(N, one(Int32))
protein_A = load_and_trans_pdb("dummy_protein.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
protein_B = load_and_trans_pdb("dummy_ligand.pdb", N)
roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
A = grid_representation(roomcoordiantes_atoms_A, N, centroids)
B = grid_representation(roomcoordiantes_atoms_B, N, centroids)
B_r = rotate_atoms(roomcoordiantes_atoms_B, rotations[1], N)
# mcA = mass_center(roomcoordiantes_atoms_A)
mcB = mass_center(B_r)
g(v::Vector3{Float32}) = (-1)*Vector3{Float32}(mcA[1]-33, mcA[2]-30, 
        mcA[3]-2) + v
h(v::Vector3{Float32}) = (-1)*Vector3{Float32}(33, 30, 2) + v
atoms_translated_A = h.(roomcoordiantes_atoms_A)
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

for i in CartesianIndices(B)
    if (B[i] != 0)
        v = Meshes.Point(i[1],i[2],i[3])
        # if (!Base.in(v, atoms_in_space_points))
            push!(atoms_in_space_points, v)
        # end
    end
end
=#
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

# show protein and rotated ligand
for i in atoms_translated_A
    v = Meshes.Point(i[1],i[2],i[3])
    push!(atoms_in_space_points, v)
end

for i in atoms_translated_B
    v = Meshes.Point(i[1],i[2],i[3])
    #if (!Base.in(v, atoms_in_space_points))
        push!(atoms_in_space_points, v)
    #end
end


viz(atoms_in_space_points, color = 1:length(atoms_in_space_points))
# viz(atoms_in_space_points)


