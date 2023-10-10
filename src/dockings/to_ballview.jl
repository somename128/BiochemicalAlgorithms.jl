using BiochemicalAlgorithms
using PDBTools

function to_ballview(pathA::String, proteinA::Molecule{Float32}, proteinB::Molecule{Float32}, rotation::Tuple, translation::Vector3{Float32}, hhb::Bool)
    
    # translate to origin and rotate + translate protein B
    atomsA = [i.r for i in eachrow(atoms_df(proteinA))]
    mcA = mass_center(atomsA) 
    trans_vecA = (-1)*Vector3{Float32}(mcA[1], mcA[2], mcA[3])
    BiochemicalAlgorithms.translate!(proteinA,trans_vecA)
    atomsB = [i.r for i in eachrow(atoms_df(proteinB))]
    mcB = mass_center(atomsB) 
    trans_vecB = (-1)*Vector3{Float32}(mcB[1], mcB[2], mcB[3])
    BiochemicalAlgorithms.translate!(proteinB,trans_vecB)
    
    # create rotation matrix for rigid body transform
    # out of rotation tuple
    rot_matrix = Matrix3(RotXYZ(rotation[1], rotation[2], rotation[3]))
    # perform translation and rotation on protein B
    rigid_transform!(proteinB, RigidTransform(rot_matrix, translation)) 

    # process to get proteinAB molecule for ballview
    sys = System()
    proteinAB = Molecule(sys, "proteinAB")  

    if (hhb) 
        atoms_A = atoms_df(proteinA)[1:1069, :]
        atoms_B = atoms_df(proteinB)[1:1123, :]
        atoms_C = atoms_df(proteinA)[1070:2138, :]
        atoms_D = atoms_df(proteinB)[1124:2246, :]

        for i in eachrow(atoms_A)
            BiochemicalAlgorithms.Atom(proteinAB, i.number, i.element, i.name, i.atom_type, i.r, i.v, i.F)
        end
        
        for i in eachrow(atoms_B)
            BiochemicalAlgorithms.Atom(proteinAB, i.number+size(atoms_A, 1)+1, i.element, i.name, i.atom_type, i.r, i.v, i.F)
        end

        for i in eachrow(atoms_C)
            BiochemicalAlgorithms.Atom(proteinAB, i.number+size(atoms_B, 1)+1, i.element, i.name, i.atom_type, i.r, i.v, i.F)
        end

        for i in eachrow(atoms_D)
            BiochemicalAlgorithms.Atom(proteinAB, i.number+size(atoms_A, 1)+size(atoms_C, 1)+2, i.element, i.name, i.atom_type, i.r, i.v, i.F)
        end
    else
        atoms_A = atoms_df(proteinA)
        atoms_B = atoms_df(proteinB)

        for i in eachrow(atoms_A)
            BiochemicalAlgorithms.Atom(proteinAB, i.number, i.element, i.name, i.atom_type, i.r, i.v, i.F)
        end
        
        for i in eachrow(atoms_B)
            BiochemicalAlgorithms.Atom(proteinAB, i.number+size(atoms_A, 1)+1, i.element, i.name, i.atom_type, i.r, i.v, i.F)
        end
    end

    # change atom coordinates in pdb of protein
    protein = readPDB(pathA)
    for i in eachindex(protein)
        protein[i].x = atoms_df(proteinAB)[i, :].r[1]
        protein[i].y = atoms_df(proteinAB)[i, :].r[2]
        protein[i].z = atoms_df(proteinAB)[i, :].r[3]
    end
    writePDB(protein,"src/dockings/ballview/2hhb/atoms_best_rotrans_grid_20.pdb")
end

hhb = true
t = Vector3{Float32}(1, 1, 1)
R = (Float32(20), Float32(20), Float32(20))
    
pathA = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
pathB = "src/dockings/testproteins/2hhb_beta_chain.pdb"
pathC = "src/dockings/ballview/2hhb.pdb"
proteinA = molecules(load_pdb(pathA))[1]
proteinB = molecules(load_pdb(pathB))[1]

to_ballview(pathC, proteinA, proteinB, R, t, hhb)