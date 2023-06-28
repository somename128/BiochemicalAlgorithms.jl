using BiochemicalAlgorithms
using FourierTools
using BenchmarkTools

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("extract_max.jl")
include("rotate_atoms.jl")
include("get_degrees.jl")

function generate_record(A::Array{Float32,3}, rotation::Matrix3{Float32}, roomcoordinates::Vector{Vector3{Float32}}, centroids::Array{Meshes.Point3, 3}, gridsize::Int32)
    # rotate atoms by rotation
    atoms = rotate_atoms(roomcoordinates, rotation, gridsize)
    # grid representation protein b
    B = grid_representation(atoms, gridsize, centroids)
    # fft-scoring
    C = ifft(fft(A).*fft(B))
    # cross-correlation over three dimensions
    # C = ccorr(A,B,3; centered=false)
    # safe α,β,γ of max fft-scoring (c)
    max = extract_max(C)
    # get degrees of rotation around x,y,z axis
    rot_in_deg = get_degrees(rotation)
    # build record for scoring table
    # transfer position of greatest value of C into
    # shifts in xyz direction
    record = (α=-((gridsize/2+max.α[1])%gridsize-gridsize/2), β=-((gridsize/2+max.β[1])%gridsize-gridsize/2), 
        γ=-((gridsize/2+max.γ[1])%gridsize-gridsize/2), R=rot_in_deg, score=max.score[1])

    return record
end