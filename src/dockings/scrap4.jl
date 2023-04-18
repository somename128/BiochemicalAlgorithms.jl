using BiochemicalAlgorithms
using JLD

include("extract_max.jl")

scoring_table = Table(α=[], β=[], γ=[], R=[], score=[])
C = load("C_matrix.jld")["C"]

max = extract_max(C)

record = (α=max.α, β=max.β, γ=max.γ, R=120, score=max.score)
push!(scoring_table,record)

# s = load("scoring_table.jld")["scoring_table"]

m = scoring_table[findmax(scoring_table.score)[2]]

