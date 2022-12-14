using BiochemicalAlgorithms
using Meshes, MeshViz
using Makie, WGLMakie

#=
point = Meshes.Point(3.5,3.5,3.5)
ball = Meshes.Ball(point, 0.5)

grid = CartesianGrid((0,0,0), (4,4,4), dims=(4,4,4))
centroids = centroid.(grid)

for i in centroids
    if(Base.in(i,ball))
        position = findall(item -> item == i, centroids)
        print(position[1])
    end
end
=#
@time begin
# load protein data from PDB
protein = load_pdb("5PTI.pdb")

# translation vector
t = Vector3{Float32}(20,20,20)

# translate protein in "positive space"
# TODO: better choice of translation vector
BiochemicalAlgorithms.translate!(protein,t)

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

# calculate mass center x,y,z coordinates
mass_x = sum(X)/length(X)
mass_y = sum(Y)/length(Y)
mass_z = sum(Z)/length(Z)

mass_center = Meshes.Point(mass_x,mass_y,mass_z)

point1 = Meshes.Point(0,0,0)
point2 = Meshes.Point(-1,0,0)

function dist(p1::Meshes.Point,p2::Meshes.Point)
    d = sqrt((p2.coords[1]-p1.coords[1])*(p2.coords[1]-p1.coords[1])+(p2.coords[2]-p1.coords[2])*(p2.coords[2]-p1.coords[2])
    +(p2.coords[3]-p1.coords[3])*(p2.coords[3]-p1.coords[3]))
    return d 
end

distances = Vector{Float64}()

 # transfer atom coordinates in mesh points
 atoms_in_space_points = Vector{Meshes.Point3}()

 for i in atoms_in_space
     v = Meshes.Point(i[1],i[2],i[3])
     push!(atoms_in_space_points, v)
 end

for i in atoms_in_space_points
    d = dist(mass_center,i)
    push!(distances, d)
end

sort(distances, rev=true)[1:100]

end