using Quaternions
using Combinatorics

include("helpers.jl")

function create_rotations2()
    
    # initialize arrays/vectors for rotation matrices and angles for rotations
    rotations = Vector{QuaternionF32}()
    # set goden ratio (1+sqrt(5)/2 ≈ 1.618)
    ϕ = Float32(1.618)
    # add all points (600) of a 4D "unit" dodecahdron (120-cell)
    # to rotations (Wikipedia: "120-cell")
    # "8"
    push!(rotations, QuaternionF32(1, 0, 0, 0))
    push!(rotations, QuaternionF32(0, 1, 0, 0))
    push!(rotations, QuaternionF32(0, 0, 1, 0))
    push!(rotations, QuaternionF32(0, 0, 0, 1))
    push!(rotations, QuaternionF32(-1, 0, 0, 0))
    push!(rotations, QuaternionF32(0, -1, 0, 0))
    push!(rotations, QuaternionF32(0, 0, -1, 0))
    push!(rotations, QuaternionF32(0, 0, 0, -1))
    # "16"
    push!(rotations, QuaternionF32(0.5, 0.5, 0.5, 0.5))
    push!(rotations, QuaternionF32(-0.5,-0.5,-0.5,-0.5))
    v = Float32[-0.5, 0.5, 0.5, 0.5]
    perm_quat!(v, rotations)
    v = Float32[-0.5, -0.5, 0.5, 0.5]
    perm_quat!(v, rotations)
    v = Float32[-0.5, -0.5, -0.5, 0.5]
    perm_quat!(v, rotations)
    # "96"
    v = Float32[0, 1/2ϕ, 0.5 , ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, -1/2ϕ, 0.5 , ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, 1/2ϕ, -0.5 , ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, -1/2ϕ, -0.5 , ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, 1/2ϕ, 0.5 , -ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, -1/2ϕ, 0.5 , -ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, 1/2ϕ, -0.5 , -ϕ/2]
    perm_quat_even!(v, rotations)
    v = Float32[0, -1/2ϕ, -0.5 , -ϕ/2]
    perm_quat_even!(v, rotations)

    return unique(rotations)
end
