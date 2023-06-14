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

protein_A = load_and_trans_pdb("src/dockings/dummy_protein.pdb", N)
roomcoordiantes_atoms_A = extract_roomcoordinates(protein_A)
protein_B = load_and_trans_pdb("src/dockings/dummy_ligand.pdb", N)
roomcoordiantes_atoms_B = extract_roomcoordinates(protein_B)
B_r = rotate_atoms(roomcoordiantes_atoms_B, rotations[4], N)
atoms_in_space_points = Base.Vector{Meshes.Point3}()
#=
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
for i in roomcoordiantes_atoms_A
    v = Meshes.Point(i[1],i[2],i[3])
    push!(atoms_in_space_points, v)
end

for i in B_r
    v = Meshes.Point(i[1],i[2],i[3])
    if (!Base.in(v, atoms_in_space_points))
        push!(atoms_in_space_points, v)
    end
end

viz(atoms_in_space_points)





