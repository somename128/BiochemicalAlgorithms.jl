using BiochemicalAlgorithms
using Meshes
using LinearAlgebra
using BenchmarkTools
using ScikitSpatial


# load protein data from PDB
protein = load_pdb("5PTI.pdb")

# translation vector
t = Vector3{Float32}(20,20,20)

# translate protein in "positive space"
# TODO: better choice of translation vector
BiochemicalAlgorithms.translate!(protein,t)

# extract room coordinates of atoms of the protein
atoms_in_space = protein.atoms.r

# radius for spheres
r = 1

# generate vector of spheres wih atom coordinates as center and radius r
spheres = Base.Vector{ScikitSpatial.Sphere}()

for a in atoms_in_space
    s = ScikitSpatial.Sphere(a,r)
    push!(spheres,s)
end

# mass center per hand
# initalize vectors for storing x y z coordinates seperately
X = Base.Vector{Float32}()
Y = Base.Vector{Float32}()
Z = Base.Vector{Float32}()

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

mass_center = Vector3{Float32}(mass_x,mass_y,mass_z)

#=
# mass_center with ScikitSpatial (not same result as mass_center per hand)
atoms_matrix = reshape(reduce(vcat, transpose.(atoms_in_space)),3,length(atoms_in_space))
mass_center_fast = ScikitSpatial.centroid(atoms_matrix)
=#

# TODO: let line rotate (three for loops might be not that clever)
surface_points = Base.Vector{Vector3{Float32}}()

for a in atoms_in_space

    # take coordinates of atom and add a bit more to find atoms "behind" every atom
    # use this point in space as endpoint of rotating line 
    # TODO: find out where point in relation to mass center is 
    endpoint = a + (100,100,100) 
    r_line = ScikitSpatial.Line(mass_center,endpoint)

    counter = 0
    max_dist = 0
    distances = Base.Vector{Tuple{Vector3{Float32},Float32}}()

    # count intersecting spheres with line, calculate distance to mass_center, safe point farest away from centroid
    # for now safe all distances too
    for i in spheres
        try
            ScikitSpatial.intersect(i, r_line)
            counter += 1
            dist = ScikitSpatial.distance(mass_center,i.point)
            push!(distances, (i.point,dist))

            if(dist > max_dist)
                max_dist = dist
            end

        catch y
            # case if line does not intersect sphere in spheres
            # throws error but should continue
            # just do nothing 
        end
    end

    # initialize vector for surface points (for only one line set at this position in code)
    # surface_points = Vector{Vector3{Float32}}()

    # distances is vector with all intersected points and their distance to mass center
    # find point which is farest away from mass center for this line (in future many lines in space) and add to surface points
    for i in distances
        if(i[2] == max_dist)
            push!(surface_points, i[1])
        end
    end

    #=
    println(counter)
    println(max_dist)
    println(distances)
    println(surface_points)
    =#
end

atoms_in_space
