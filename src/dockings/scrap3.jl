using BiochemicalAlgorithms
using BenchmarkTools
using Makie, WGLMakie
using Meshes, MeshViz
using JLD
using ProfileView
using ProgressBars
using Profile
using TypedTables
using LinearAlgebra
using FFTW
using FourierTools
using Rotations

include("grid_representation.jl")
include("load_trans_pdb.jl")
include("min_max_atoms.jl")
include("create_rotations.jl")
include("set_marked_cells.jl")
include("generate_record.jl")
include("extract_roomcoordinates.jl")
include("create_atomballs.jl")
include("correlation_docking.jl")
include("mass_center.jl")
include("rotate_atoms.jl")
include("helpers.jl")
include("set_gridsize.jl")
include("helpers.jl")

N = set_gridsize("dummy_protein_vol2.pdb","dummy_ligand_vol2.pdb")
a = Int32(31)
interp!(a,N)
