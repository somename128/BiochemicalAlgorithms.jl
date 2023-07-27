using PDBTools
# script for creating a 10x10x10 cube with a 5x5x5 cube missing in
# the lower left corner (protein) and a 6x6x6 cube (ligand)
# for the purpose of easy geometries 
# all the stuff for still using pdb instead of only an array 

# initialize all sets and arrays (sets because difference between two sets needed - protein)
atoms_protein = Set()
atoms_ligand = []
atoms_remove = Set()
atoms_protein_cleaned = Set()
atoms_protein_indexed = []
# sizes of the cubes
NP = 10
NL = 6
NR = 5

# fill sets/arrays with points inside hexagons of grid
# to have no overlap and exactly the amount of colored cells 
# during algorithm
# portein cube
for i in 1:NP+1, j in 1:NP+1, k in 1:NP+1
  v = [i,j,k]
  push!(atoms_protein, v)
end

# ligand cube
for i in 1:NL+1, j in 1:NL+1, k in 1:NL+1
  v = [i,j,k]
  push!(atoms_ligand, v)
end

# missing part in protein cube
for i in 1:NR+1, j in 1:NR+1, k in 1:NR+1
  v = [i,j,k]
  push!(atoms_remove, v)
end

# create "missing edge" in portein cube 
atoms_protein_cleaned = setdiff(atoms_protein, atoms_remove)

# sort the points 
sort!(collect(atoms_protein_cleaned))

# transfer set to array/vector
for i in atoms_protein_cleaned
  push!(atoms_protein_indexed,i)
end


# change atom coordinates in pdb of protein
protein = readPDB("2ptc_protein_copy.pdb")
for i in eachindex(atoms_protein_indexed)
   protein[i].x = atoms_protein_indexed[i][1]
   protein[i].y = atoms_protein_indexed[i][2]
   protein[i].z = atoms_protein_indexed[i][3]
   protein[i].segname = string(protein[i].name[1])
end

# change atom coordinates in pdb of ligand
ligand = readPDB("2ptc_ligand_copy.pdb")
for i in eachindex(atoms_ligand)
  ligand[i].x = atoms_ligand[i][1]
  ligand[i].y = atoms_ligand[i][2]
  ligand[i].z = atoms_ligand[i][3]
  ligand[i].segname = string(ligand[i].name[1])
end

# write changed atoms of protein and/or ligand (only atoms of pdb) in pdb file
writePDB(protein,"atoms_dummy_protein_vol4.pdb")
writePDB(ligand,"atoms_dummy_ligand_vol4.pdb")

  