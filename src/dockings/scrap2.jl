using Distributions
using Rotations
using Quaternions
using BenchmarkTools
using JLD2
using SpecialFunctions

include("bingham_functions.jl")
include("create_rotations.jl")
include("correlation_docking.jl")

results_docking = load_object("results_docking.jld2")

scoring_table = results_docking[1]
# extract N best quaternions (current value: twenty)
Q = Vector{QuaternionF32}()
for i in eachrow(scoring_table)
    q = extract_quaternion(i)
    push!(Q, q)
end

# generate lookup table F for Bingham dirstribution
F = create_lookup_table_F()
# estimate V
# V = estimate_V(Q)
# estimate Λ
# Λ = estimate_Λ(Q, V, F)

#=
# generate scatter matrix
S = scatter_matrix(Q)
# sample new point
R = metropolis_hastings_sampler(Q[1], Λ, V, S, F, Int32(10))

R *= 1/norm(R)
norm(R)
=#

