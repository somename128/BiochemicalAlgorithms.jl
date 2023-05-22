using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using ProgressBars
using Profile
using TypedTables

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")
include("generate_record.jl")
include("reset_rotation.jl")

N = 128
protein = load_and_trans_pdb("2ptc_protein.pdb",N)
atoms_before = Vector{Vector{Float32}}()
atoms_after = Vector{Vector{Float32}}()
atoms_transpose = Vector{Vector{Float32}}()
# centroids = create_centroids(N,1)
# A = grid_representation(protein,N, centroids)
rotations = create_rotations()
for i in atoms_df(protein).r
    push!(atoms_before,i)
end
rigid_transform!(protein, rotations[5])
for i in atoms_df(protein).r
    push!(atoms_after,i)
end
reset_rotation!(protein,rotations[5])
for i in atoms_df(protein).r
    push!(atoms_transpose,i)
end
println(atoms_before ≈ atoms_after)
println(atoms_before ≈ atoms_transpose)
#=
t = Vector3{Float32}(0,0,0)
R = Matrix3{Float32}([0 0 0; 
            0 0 0; 
            0 0 0])
rt = RigidTransform{Float32}(R,t)


table = Table(α=[0.0], β=[0.0], γ=[0.0], R=[rt], score=[0.0])


# Profile.clear()
@time record = generate_record(A, rotations[5], "2ptc_ligand.pdb", centroids,N)
# ProfileView.view()

table[1] = (α=record.α, β=record.β, γ=record.γ, R=record.R, score=record.score)
=#
