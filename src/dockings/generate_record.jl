using FFTW

include("grid_representation.jl")
include("extract_max.jl")
include("rotate_atoms.jl")
include("get_degrees.jl")
include("extract_max_fast.jl")

function generate_record(A::Array{ComplexF32,3}, rotation::QuaternionF32, roomcoordinates::Vector{Tuple{String, Vector3{Float32}}}, centroids::Array{Meshes.Point3f, 3}, gridsize::Int32, res::Int32, vdW::Bool)
    # rotate atoms roomcoordinates by rotation
    atoms = rotate_atoms(roomcoordinates, rotation, gridsize)
    # grid representation protein b
    B = grid_representation(atoms, gridsize, centroids, res, true, vdW)
    # fft-scoring
    C = ifft(fft(A).*conj(fft(B)))
    # safe α,β,γ of max fft-scoring (c)
    max = extract_max_fast(C)
    
    # build record for scoring table
    # transfer position of greatest value of C into shifts 
    # in xyz direction
    # minus 1/res in each direction for start at (0,0,0) in upper 
    # left corner (not (1/res,1/res,1/res)) (this calculation is
    # in helper function in helper.jl)

    return (α=Float32(interp!(max[1],gridsize,res)), β=Float32(interp!(max[2],gridsize,res)), 
    γ=Float32(interp!(max[3],gridsize,res)), R=get_degrees(rotation), score=Float32(max[4]))
end