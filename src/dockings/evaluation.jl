using BiochemicalAlgorithms
using Rotations

function evaluation(protein_A::String, protein_B::String, protein_complex::String, rotation::Tuple, translation::Vector3{Float32})
    # load protein A,B and complex AB
    proteinA = molecules(load_pdb(protein_A))[1]
    proteinB = molecules(load_pdb(protein_B))[1]
    complexAB = molecules(load_pdb(protein_complex))[1]

    # create rotation matrix for rigid body transform
    # out of rotation tuple
    rot_matrix = Matrix3(RotXYZ(rotation[1], rotation[2], rotation[3]))
    # create rigid transform for function 
    rot_and_trans = RigidTransform(rot_matrix, translation)
    # perform translation and rotation on protein B
    rigid_transform!(proteinB, rot_and_trans)    

    # process to get proteinAB molecule and docked protein molecule 
    # for rmsd
    proteinAB = Molecule("proteinAB")
    atoms_AB = unique(atoms_df(complexAB), 2)
    docked_proteinAB = Molecule("docked_proteinAB")    
    atoms_A = unique(atoms_df(proteinA), 2)
    atoms_B = unique(atoms_df(proteinB), 2)
    
    for i in eachrow(atoms_AB)
        Atom(proteinAB, i.number, i.element, i.name, i. atom_type, i.r, i.v, i.F)
    end

    for i in eachrow(atoms_A)
        Atom(docked_proteinAB, i.number, i.element, i.name, i. atom_type, i.r, i.v, i.F)
    end
    
    for i in eachrow(atoms_B)
        Atom(docked_proteinAB, i.number, i.element, i.name, i. atom_type, i.r, i.v, i.F)
    end
    
    return compute_rmsd(TrivialAtomBijection(docked_proteinAB, proteinAB))

end