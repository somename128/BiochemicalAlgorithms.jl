using PDBTools


atoms_protein = Set()
atoms_ligand = []
atoms_remove = Set()
atoms_protein_cleaned = Set()
atoms_protein_indexed = []
NP = 10
NL = 6
NR = 5

for i in 2:NP+2, j in 2:NP+2, k in 2:NP+2
  v = [i+0.5,j+0.5,k+0.5]
  push!(atoms_protein, v)
end

for i in 2:NL+2, j in 2:NL+2, k in 2:NL+2
  v = [i+0.5,j+0.5,k+0.5]
  push!(atoms_ligand, v)
end

for i in 2:NR+2, j in 2:NR+2, k in 2:NR+2
  v = [i+0.5,j+0.5,k+0.5]
  push!(atoms_remove, v)
end

atoms_protein_cleaned = setdiff(atoms_protein, atoms_remove)

sort!(collect(atoms_protein_cleaned))

for i in atoms_protein_cleaned
  push!(atoms_protein_indexed,i)
end

protein = readPDB("2ptc_protein_copy.pdb")


for i in eachindex(atoms_protein_indexed)
   protein[i].x = atoms_protein_indexed[i][1]
   protein[i].y = atoms_protein_indexed[i][2]
   protein[i].z = atoms_protein_indexed[i][3]
   protein[i].segname = string(protein[i].name[1])
end

writePDB(protein,"atoms_dummy_protein.pdb")

  