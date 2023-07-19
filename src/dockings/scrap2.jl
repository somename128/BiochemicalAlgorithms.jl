using Distributions
using Rotations
using Quaternions

include("quaternion_functions.jl")

q1 = quat_from_axisangle([1,0,0], deg2rad(20))

