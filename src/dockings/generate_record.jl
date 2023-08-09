using FFTW

include("grid_representation.jl")
include("extract_max.jl")
include("rotate_atoms.jl")
include("get_degrees.jl")

function generate_record(A::Array{Float32,3}, rotation::QuaternionF32, roomcoordinates::Vector{Tuple{String, Vector3{Float32}}}, centroids::Array{Meshes.Point3f, 3}, gridsize::Int32, res::Int32)
    # rotate atoms roomcoordinates by rotation
    atoms = rotate_atoms(roomcoordinates, rotation, gridsize)
    # grid representation protein b
    B = grid_representation(atoms, gridsize, centroids, res, true, false)
    # fft-scoring
    C = ifft(fft(A).*conj(fft(B)))
    # safe α,β,γ of max fft-scoring (c)
    max = extract_max(C)
    # get degrees of rotation around x,y,z axis
    rot_in_deg = get_degrees(rotation)
    # build record for scoring table
    # transfer position of greatest value of C into shifts 
    # in xyz direction
    # minus 1 in each direction for start at (0,0,0) in upper 
    # left corner (not (1,1,1))
    record = (α=interp!(max.α[1],gridsize,res)-one(Float32), β=interp!(max.β[1],gridsize,res)-one(Float32), 
        γ=interp!(max.γ[1],gridsize,res)-one(Float32), R=rot_in_deg, score=max.score[1])

    return record
end