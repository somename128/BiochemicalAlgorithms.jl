using Meshes, MeshViz
using Makie, WGLMakie

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