using Quaternions

include("helpers.jl")

function create_rotations2()
    
    # initialize arrays/vectors for rotation matrices and angles for rotations
    rotations = Vector{QuaternionF32}()
    # set goden ratio (1+sqrt(5)/2 ≈ 1.618)
    ϕ = Float32(1.618)
    # calculate uniformly distributed points on a unit 4D hypersphere
    # by collecting/calculating all 600 vertices of a hyperdodecahedron
    # (120-cell)
    # Paper: Mamone, Levitt "Orientational Sampling Schemes Based on Four
    # Dimensional Polytopes (p. 1442 Table 3)
    # collect all points of a hypertetrahedron (5-cell)
    S = Vector{QuaternionF32}()
    push!(S, QuaternionF32(1, 0, 0, 0))
    push!(S, QuaternionF32(-1/4, sqrt(5)/4, sqrt(5)/4, sqrt(5)/4))
    push!(S, QuaternionF32(-1/4, -sqrt(5)/4, -sqrt(5)/4, sqrt(5)/4))
    push!(S, QuaternionF32(-1/4, -sqrt(5)/4, sqrt(5)/4, -sqrt(5)/4))
    push!(S, QuaternionF32(-1/4, sqrt(5)/4, -sqrt(5)/4, -sqrt(5)/4))
    # collect all points of a hypericosahedron (600-cell)
    I = Vector{QuaternionF32}()
    # 8-cell
    v = Float32[1, 0, 0, 0]
    perm_quat_even!(v, I)
    v = Float32[-1, 0, 0, 0]
    perm_quat_even!(v, I)
    # 16-cell
    v = Float32[0.5, 0.5, 0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, 0.5, 0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, 0.5, -0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, 0.5, -0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, -0.5, 0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, -0.5, 0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, -0.5, -0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[0.5, -0.5, -0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, 0.5, 0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, 0.5, 0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, -0.5, 0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, 0.5, -0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, -0.5, 0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, -0.5, 0.5, -0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, -0.5, -0.5, 0.5]
    perm_quat_even!(v, I)
    v = Float32[-0.5, -0.5, -0.5, -0.5]
    perm_quat_even!(v, I)
    # remaining points
    v = Float32[ϕ/2, 0.5, 1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[ϕ/2, 0.5, -1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[ϕ/2, -0.5, 1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[ϕ/2, -0.5, -1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[-ϕ/2, 0.5, 1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[-ϕ/2, 0.5, -1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[-ϕ/2, -0.5, 1/2ϕ, 0]
    perm_quat_even!(v, I)
    v = Float32[-ϕ/2, -0.5, -1/2ϕ, 0]
    perm_quat_even!(v, I)

    # calculate all 600 points of 120-cell
    # description in paper
    # multply q0 with all possible quaternion products of S and I
    # multiple same vectors in I therefore unique()
    q0 = QuaternionF32(1/sqrt(2), 1/sqrt(2), 0, 0)
    for q1 in S, q2 in unique(I)
        push!(rotations, q0*q1*q2)
    end

    return rotations
end
