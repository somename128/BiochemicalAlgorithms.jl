# functions from quaternions.jl
using Quaternions
using LinearAlgebra

# build the quaternion which represents rotation around a axis 
function quat_from_axisangle(axis::AbstractVector, theta::Real)
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    s, c = sincos(theta / 2)
    axis = normalize(axis)
    return QuaternionF32(c, s*axis[1], s*axis[2], s*axis[3])
end

# recovers the axis and the angle from a quaternion
function axisangle_from_quat(q::QuaternionF32)
    axis = Vector3{Float32}(q.v1, q.v2, q.v3)/sqrt(q.v1^2 + q.v2^2 + q.v3^2)
    angle = 2*atan(sqrt(q.v1^2 + q.v2^2 + q.v3^2), q.s)
    
    return axis, angle 
end

# rotate vector with quaternion
function rotate_vector(q::QuaternionF32, u::AbstractVector)
    if length(u) != 3
        error("Must be a 3-vector")
    end
    q_u = QuaternionF32(0, u[1], u[2], u[3])
    q_v = q*q_u*conj(q)
    return Vector3{Float32}(imag_part(q_v))
end

# convert quaternion to rotation matrix
function rotmatrix_from_quat(q::QuaternionF32)
    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3
    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3
    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2
    r = [1 - (yy + zz)     xy - sz     xz + sy;
            xy + sz   1 - (xx + zz)    yz - sx;
            xz - sy      yz + sx  1 - (xx + yy)]
    return r
end

# convert rotation matrix to quaternion
function quat_from_rotmatrix(dcm::AbstractMatrix{T}) where {T<:Real}
    a2 = 1 + dcm[1,1] + dcm[2,2] + dcm[3,3]
    a = sqrt(a2)/2
    b,c,d = (dcm[3,2]-dcm[2,3])/4a, (dcm[1,3]-dcm[3,1])/4a, (dcm[2,1]-dcm[1,2])/4a
    return QuaternionF32(a,b,c,d)
end