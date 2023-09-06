using Meshes

p = Meshes.Point3(0.5,0.5,0.5)

b = Meshes.Ball((0,0,0), 0.75)

Base.in(p,b)
