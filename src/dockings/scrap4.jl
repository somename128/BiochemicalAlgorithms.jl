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
atoms_in_space_points = Base.Vector{Meshes.Point3}()

for i in eachindex(rotations)
    B_r = rotate_atoms(roomcoordiantes_atoms_B,rotations[i])

    for a in B_r
        v = Meshes.Point(a[1],a[2],a[3])
        if (!Base.in(v, atoms_in_space_points))
            push!(atoms_in_space_points, v)
        end
    end
    println(i,"/",length(rotations))
end


#=
centroids = create_centroids(N, one(Int32))
A = grid_representation(roomcoordiantes_atoms_A, N, centroids)
B = grid_representation(roomcoordiantes_atoms_B, N, centroids)
C = ifft(fft(A).*fft(B))
C_real = real(C)
C_sparse = zero_small!(C_real,0.5)

atoms_in_space_points = Base.Vector{Meshes.Point3}()
for i in roomcoordiantes_atoms_A
    v = Meshes.Point(i[1],i[2],i[3])
    push!(atoms_in_space_points, v)
end

for i in B_r1
    v = Meshes.Point(i[1],i[2],i[3])
    if (!Base.in(v, atoms_in_space_points))
        push!(atoms_in_space_points, v)
    end
end
=#
viz(atoms_in_space_points)



 
#=
min_max = min_max_atoms(atoms)
atoms_r = rotate_atoms(atoms,rotations[2034])
min_max_atoms(atoms_r)
centroids = create_centroids(N,one(Int32))
atomballs = create_atomballs(atoms_r)
colored_cells = set_marked_cells(atomballs, centroids, atoms_r)
grid = create_inner_outer_grid(colored_cells, N)
=#
# transfer atom coordinates in mesh points


#=
for i in CartesianIndices(grid)
if(grid[i] != 0)
    v = Meshes.Point(i[1],i[2],i[3])
    push!(atoms_in_space_points, v)
end
end
=#




