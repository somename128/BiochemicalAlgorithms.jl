using BiochemicalAlgorithms
using Meshes
using Quaternions
using BenchmarkTools
using Combinatorics


include("load_trans_pdb.jl")
include("extract_roomcoordinates.jl")
include("rotate_atoms.jl")
include("create_rotations.jl")
include("bingham_functions_vol2.jl")
include("helpers.jl")
include("create_rotations2.jl")

rotations = create_rotations2()


