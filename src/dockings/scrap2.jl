using DelimitedFiles
using Makie, WGLMakie
using Meshes, MeshViz

include("load_trans_pdb.jl")
include("set_gridsize.jl")
include("grid_representation.jl")
include("extract_roomcoordinates.jl")
include("create_centroids.jl")

path = "simple_geometry/cube_origin_vdW.pdb"
res = Int32(2)
N = set_gridsize(path, path)
protein_origin = molecules(load_pdb(path))[1]
protein = load_and_trans_pdb(path, N)
centroids = create_centroids(N, res)
coord = extract_roomcoordinates(protein)
cube = grid_representation(coord, N, centroids, res, false, true)
# write txt file
#=
f = open("cube_origin_vdW.txt", "a")
for k in 4:13
    write(f, string(k,".layer\n"))
    flush(f)
    for i in 3:14
        for j in 3:14
            write(f, string(cube[i,j,k],"  "))
            flush(f)   
        end
    write(f, "\n")
    flush(f)
    end
end

close(f)
=#

# visualize points 
coord_origin = extract_roomcoordinates(protein_origin)
atoms = [convert(Meshes.Point3f, i[2]) for i in coord]
# cents = Vector{Meshes.Point3f}()
# for i in centroids
#      push!(cents, i)
# end
#=
atoms = Base.Vector{Meshes.Point3f}()
for i in CartesianIndices(cube)
    if (cube[i] != 0 && cube[i] != -15)
        v = Meshes.Point(i[1]/res-1/2res,i[2]/res-1/2res,i[3]/res-1/2res)
        push!(atoms, v)
    end
end
=#
# without color
# viz(atoms)
# with color
viz(atoms, color = 1:length(atoms))
