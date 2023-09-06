using PDBTools
# script for creating a 10x10x10 cube with a 5x5x5 cube missing in
# the lower left corner (protein) and a 6x6x6 cube (ligand)
# for the purpose of easy geometries 
# all the stuff for still using pdb instead of only an array 

# initialize all sets and arrays (sets because difference between two sets needed - protein)
atoms_cube = []

# fill sets/arrays with points inside hexagons of grid
# to have no overlap and exactly the amount of colored cells 
# during algorithm
# portein cube
for i in 1:10, j in 1:10, k in 1:10
    v = [(i-1)-4.5,(j-1)-4.5,(k-1)-4.5]
    push!(atoms_cube, v)
end


# change atom coordinates in pdb of ligand
ligand = readPDB("2ptc_protein_copy.pdb")
for i in eachindex(atoms_cube)
  ligand[i].x = atoms_cube[i][1]
  ligand[i].y = atoms_cube[i][2]
  ligand[i].z = atoms_cube[i][3]
  ligand[i].segname = "D"
end

# write changed atoms of protein and/or ligand (only atoms of pdb) in pdb file
writePDB(ligand,"atoms_cube_origin_huge.pdb")

