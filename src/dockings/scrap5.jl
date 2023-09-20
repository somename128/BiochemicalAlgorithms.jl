using Meshes, MeshViz
using Makie, WGLMakie
using ProfileView
using Profile
using BenchmarkTools
using JET
using BiochemicalAlgorithms
using FFTW

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("create_rotations.jl")
include("create_centroids.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("set_gridsize.jl")
include("create_rotations2.jl")
include("create_atomballs.jl")
include("set_surface_cells_fast.jl")
include("set_surface_cells.jl")
include("set_marked_cells_fast.jl")
include("set_marked_cells.jl")
include("rotate_atoms.jl")
include("create_inner_outer_grid.jl")
include("transform_cells!.jl")
include("set_grid.jl")
include("set_grid_vdW.jl")


pathA = "src/dockings/testproteins/2hhb_beta_chain.pdb"
pathB = "src/dockings/testproteins/2ptc_inhibitor.pdb"
# init
res = Int32(1)
# N = set_gridsize(pathA, pathB)
N = Int32(150)
println("Gridsize: ", N)
# N = Int32(32)
# load and translate protein a
protein_A = load_and_trans_pdb(pathA, N)
atoms = extract_roomcoordinates(protein_A)
# calculate centroids in a NxNxN grid with cells
# of 1 angström, only done once
centroids = create_centroids(N, res)
thickness = Float32(2)
inner_radius = create_atomballs(atoms, -(thickness/2))
outer_radius = create_atomballs(atoms, thickness/2)
roomcoordinates = Array{Vector3{Float32}}(undef, length(atoms))
[roomcoordinates[i] = atoms[i][2] for i in eachindex(atoms)]
rotations = create_rotations()

is_smaller = true
inner = set_marked_cells(inner_radius, centroids, roomcoordinates, res)
surf = set_surface_cells(inner_radius, outer_radius, centroids, roomcoordinates, res)
grid1 = create_inner_outer_grid_vdW(inner, surf, N, res, is_smaller)
grid2 = set_grid_vdW(inner_radius, outer_radius, centroids, roomcoordinates, N, res, is_smaller)

@time begin
inner = set_marked_cells(inner_radius, centroids, roomcoordinates, res)
surf = set_surface_cells(inner_radius, outer_radius, centroids, roomcoordinates, res)
grid1 = create_inner_outer_grid_vdW(inner, surf, N, res, is_smaller)
end
@time grid2 = set_grid_vdW(inner_radius, outer_radius, centroids, roomcoordinates, N, res, is_smaller)
count = 0
count1 = 0
count2 = 0
for i in eachindex(grid1)
    if (grid1[i] != grid2[i])
        println("position: ", i, ", value old: ", grid1[i], ", value new: ", grid2[i])
        global count += 1
    end

    if (grid1[i] != 0)
        global count1 += 1
    end

    if (grid2[i] != 0)
        global count2 += 1
    end

end 
println(count1)
println(count2)
count