using Meshes
using LinearAlgebra
using ScikitSpatial


s1 = ScikitSpatial.Sphere([0,0,0],1)
s2 = ScikitSpatial.Sphere([10,20,30],1)
s3 = ScikitSpatial.Sphere([0,0,0],1)
spheres = [s1,s2,s3]
line = ScikitSpatial.Line([0,0,0],[1,1,1])

intersecting_spheres = 0

for s in spheres
    try
        ScikitSpatial.intersect(s,line)

        if(!isempty(ScikitSpatial.intersect(s,line)))
            global intersecting_spheres += 1
        end
    catch y
    end
end

print(intersecting_spheres)








