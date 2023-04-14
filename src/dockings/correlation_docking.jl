using BiochemicalAlgorithms
using BenchmarkTools
using FFTW
using JLD

include("load_trans_pdb.jl")
include("grid_representation.jl")

@time begin
# load and translate protein a
protein_A = load_and_trans_pdb("2ptc_protein.pdb")
# load and translate protein b
protein_B = load_and_trans_pdb("2ptc_ligand.pdb")

# grid representation protein a
A = grid_representation(protein_A)

# rotate protein b by R
# grid representation protein b
B = grid_representation(protein_B)


# fft-scoring
C = ifft(fft(A).*fft(B))
end

@time begin
save("C_matrix.jld", "C", C)
end
# safe α,β,γ,R of max fft-scoring (c)

# loop to rotation if "greatest" c not reached
