using BiochemicalAlgorithms

protein = load_pdb("2ptc_protein.pdb")

# translate protein in positive space
println("Translate protein in positive space...")
# extract room coordinates of atoms of the protein
atoms_in_space = protein.atoms.r
# initalize vectors for storing x y z coordinates seperately
X = Vector{Float32}()
Y = Vector{Float32}()
Z = Vector{Float32}()

# fill vectors with coordinates
for i in atoms_in_space
    push!(X, i[1])
    push!(Y, i[2])
    push!(Z, i[3])
end

# calculate x y z min for translation (max maybe needed one time)
min_x = minimum(X)
max_x = maximum(X)
min_y = minimum(Y)
max_y = maximum(Y)
min_z = minimum(Z)
max_z = maximum(Z)

DIA = Vector{Float32}()
push!(DIA,abs(min_x-max_x)/2)
push!(DIA,abs(min_y-max_y)/2)
push!(DIA,abs(min_z-max_z)/2)

maximum(DIA)
