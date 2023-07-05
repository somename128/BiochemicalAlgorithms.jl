using Distributions
using Rotations

include("create_rotations.jl")

rotations = create_rotations()

# get random new numbers for rotations
α = deg2rad(rand(Uniform(Float32(-20),Float32(20))))
β = deg2rad(rand(Uniform(Float32(-20),Float32(20))))
γ = deg2rad(rand(Uniform(Float32(-20),Float32(20))))

# create rotation and push in vector
r = RotXYZ{Float32}(α,β,γ)

v = Vector3{Float32}(1,2,3)

v_rotated = r*v

println(r)
r[1,3]
