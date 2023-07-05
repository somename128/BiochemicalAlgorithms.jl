include("correlation_docking.jl")

@time score = correlation_docking("src/dockings/dummy_protein_vol2.pdb","src/dockings/dummy_ligand_vol3.pdb")[1]

