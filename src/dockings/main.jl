include("correlation_docking.jl")

@time score = correlation_docking("src/dockings/dummy_protein.pdb","src/dockings/dummy_ligand.pdb")[1]

