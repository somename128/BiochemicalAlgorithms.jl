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
for i in 1:NP, j in 1:NP, k in 1:NP
  v = [(i-1)-4.5,(j-1)-4.5,(k-1)-4.5]
  push!(atoms_protein, v)
end

# ligand cube
for i in 1:NL, j in 1:NL, k in 1:NL
  v = [(i-1)-2.5,(j-1)-2.5,(k-1)-2.5]
  push!(atoms_ligand, v)
end

# missing part in protein cube
for i in 1:NR, j in 1:NR, k in 1:NR
  v = [(i-1)-4.5,(j-1)-4.5,(k-1)-4.5]
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
   protein[i].segname = "D"
end

# change atom coordinates in pdb of ligand
ligand = readPDB("2ptc_ligand_copy.pdb")
for i in eachindex(atoms_ligand)
  ligand[i].x = atoms_ligand[i][1]
  ligand[i].y = atoms_ligand[i][2]
  ligand[i].z = atoms_ligand[i][3]
  ligand[i].segname = "D"
end

# write changed atoms of protein and/or ligand (only atoms of pdb) in pdb file
writePDB(protein,"atoms_cube_origin_huge_A.pdb")
writePDB(ligand,"atoms_cube_origin_huge_B.pdb")

  