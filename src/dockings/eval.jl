using BiochemicalAlgorithms
using Rotations

include("mass_center.jl")
include("helpers.jl")

function eval(protein_A::Molecule{Float32}, protein_B::Molecule{Float32}, protein_complex::Molecule{Float32}, rotation::Tuple, translation::Vector3{Float32})
 
    # translate protein A and B in origin
    # see load_and_trans_pdb.jl
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
    # create rigid transform for function 
    rot_and_trans = RigidTransform(rot_matrix, translation)
    # perform translation and rotation on protein B
    rigid_transform!(proteinB, rot_and_trans)    

    # process to get proteinAB molecule for rmsd
    sys = System()
    proteinAB = Molecule(sys, "proteinAB")   
    atoms_A = atoms_df(proteinA)
    atoms_B = atoms_df(proteinB)

    for i in eachrow(atoms_A)
        BiochemicalAlgorithms.Atom(proteinAB, i.number, i.element, i.name, i.atom_type, i.r, i.v, i.F)
    end
    
    for i in eachrow(atoms_B)
        BiochemicalAlgorithms.Atom(proteinAB, i.number+size(atoms_A, 1)+1, i.element, i.name, i.atom_type, i.r, i.v, i.F)
    end

    # map_rigid! makes best rotation for mapping to compute rmsd
    map_rigid!(proteinAB, complexAB)

    return compute_rmsd(TrivialAtomBijection(proteinAB, complexAB))

end