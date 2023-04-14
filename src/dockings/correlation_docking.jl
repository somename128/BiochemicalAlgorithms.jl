using BiochemicalAlgorithms
using BenchmarkTools
using FFTW
using JLD
using TypedTables

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("helpers.jl")

#=
# load and translate protein a
protein_A = load_and_trans_pdb("2ptc_protein.pdb")
# load and translate protein b
protein_B = load_and_trans_pdb("2ptc_ligand.pdb")

# grid representation protein a
# A = grid_representation(protein_A)

# rotate protein b by R
# grid representation protein b
B = grid_representation(protein_B)

# fft-scoring
C = ifft(fft(A).*fft(B))

save("C_matrix.jld", "C", C)
=#
C = real(load("C_matrix.jld")["C"])
# safe α,β,γ,R of max fft-scoring (c)
t = Table(x=[], y= [], z=[], score=[])

C_sparse = zero_small!(C,0.5)

for k = 1:size(C_sparse,3)
    for j = 1:size(C_sparse,2)
        for i = 1:size(C_sparse,1)
            push!(t,(x=i, y=j, z=k, score=C_sparse[i,j,k]))
        end
    end
end

t[findmax(real(t.score))[2]]

# loop to rotation if "greatest" c not reached

