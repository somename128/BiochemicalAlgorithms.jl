using BiochemicalAlgorithms
using BenchmarkTools
using FFTW
using JLD
using TypedTables

include("load_trans_pdb.jl")
include("grid_representation.jl")
include("extract_max.jl")
include("create_rotations.jl")

scoring_table = Table(α=[], β= [], γ=[], score=[])

# load and translate protein a
# protein_A = load_and_trans_pdb("2ptc_protein.pdb")
# load and translate protein b
# protein_B = load_and_trans_pdb("2ptc_ligand.pdb")
# grid representation protein a
# A = grid_representation(protein_A)
# save("grid_repo_A.jld", "A", A)
A = load("grid_rep_A.jld")["A"]

# get rotations - build via rigid_transform! with translation vector (0,0,0)
rotations = create_rotations()
# rotate protein b by R
for R in rotations
    # load and translate protein b per loop to use rotations relative to 
    # origin coordinates
    protein_B = load_and_trans_pdb("2ptc_ligand.pdb")
    rigid_transform!(protein_B, R)
    # grid representation protein b
    B = grid_representation(protein_B)
    # fft-scoring
    C = ifft(fft(A).*fft(B))
    # -----------------------------------------
    # only for time saving purposes, not part of the real algorithm
    # save("C_matrix.jld", "C", C)

    # C = load("C_matrix.jld")["C"]
    # -----------------------------------------
    # safe α,β,γ,R of max fft-scoring (c)
    max = extract_max(C)
    push!(scoring_table,max)
end    
# --------------------------------------------
# safe results of run
save("scoring_table.jld", "scoring_table", scoring_table)
# --------------------------------------------
scoring_table

# loop to rotation if "greatest" c not reached

