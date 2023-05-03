using BiochemicalAlgorithms
using Makie, WGLMakie
using Meshes, MeshViz
using LinearAlgebra
using BenchmarkTools

include("load_trans_pdb.jl")
include("helpers.jl")
include("min_max_atoms.jl")
include("mass_center.jl")
include("create_centroids.jl")

sys = load_pdb("2ptc_ligand.pdb")

mol = molecules(sys)[1]

println(atoms_df(mol).r[1])

println(mass_center(mol))



protein = load_and_trans_pdb("2ptc_ligand.pdb",128)

t = Vector3{Float32}(0,0,0)

Θ = 5
ψ = 5
ϕ = 5

R = Matrix3{Float32}([cosd(Θ)*cosd(ϕ) cosd(Θ)*sind(ϕ) -sind(Θ); 
            sind(ψ)*sind(Θ)*cosd(ϕ)-cosd(ψ)*sind(ϕ) sind(ψ)*sind(Θ)*sind(ϕ)+cosd(ψ)*cosd(ϕ) cosd(Θ)*sind(ψ); 
            cosd(ψ)*sind(Θ)*cos(ϕ)+sind(ψ)*sind(Θ) cosd(ψ)*sind(Θ)*sind(ϕ)-sind(ψ)*cosd(ϕ) cosd(Θ)*cosd(ψ)])

rigidtransform = RigidTransform(R,t)

@time rigid_transform!(protein,rigidtransform)





