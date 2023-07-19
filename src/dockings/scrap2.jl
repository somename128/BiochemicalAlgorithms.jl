using Distributions
using Rotations
using Quaternions

include("create_gridpoints.jl")

gridpoints = create_gridpoints(Int32(32), Int32(4))