using BiochemicalAlgorithms

include("load_trans_pdb.jl")

protein_B = load_and_trans_pdb("2ptc_ligand.pdb")

for i in 1:2
    protein_B = load_and_trans_pdb("2ptc_ligand.pdb")
    display(protein_B)

    t = Vector3{Float32}(0,0,0)
    R = Matrix3{Float32}([10 10 10; 10 10 10; 10 10 10])

    rotation = RigidTransform(R,t)

    rigid_transform!(protein_B, rotation)
end