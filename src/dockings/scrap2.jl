using DelimitedFiles
using Makie, WGLMakie
using Meshes, MeshViz

include("load_trans_pdb.jl")
include("set_gridsize.jl")
include("grid_representation.jl")
include("extract_roomcoordinates.jl")
include("create_centroids.jl")

path = "src/dockings/simple_geometry/cube_origin.pdb"
res = Int32(1)
N = set_gridsize(path, path)
println(N)
protein_origin = molecules(load_pdb(path))[1]
protein = load_and_trans_pdb(path, N)
centroids = create_centroids(N, res)
coord = extract_roomcoordinates(protein)
cube = grid_representation(coord, N, centroids, res, false, false)
# write txt file
#=
f = open("dot_origin_vdW_2.txt", "a")
for k in 1:N*res
    write(f, string(k,".layer\n"))
    flush(f)
    for i in 1:N*res
        for j in 1:N*res
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
atoms = [convert(Meshes.Point3f, i[2]) for i in coord_origin]
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
