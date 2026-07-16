using FerriteGmsh

gmsh.initialize()
gmsh.model.add("battery_3d")

# definitions
r1 = .1
r2 = .125
h1 = .2
h2 = .25
ri = .04

d = (h2 - h1)/2

# define namespace for easier typing
occ = gmsh.model.occ

# cylinders
cyl1 = occ.addCylinder(0, 0, d, 0, 0, h1, r1)
cyl2 = occ.addCylinder(0, 0, 0, 0, 0, d, ri)

# merge
inner = occ.fuse([(3, cyl1)], [(3, cyl2)])
outer = occ.addCylinder(0, 0, 0, 0, 0, h2, r2)

frag = factory.fragment([(3, outer)], [inner[1][1]])

# frag[1] gives new entities from the fragment function
# frag[1][1] gives first (dim, tag)-tuple of the old fragment function
# frag[1][1][2] gives tag of outer surface
gmsh.model.addPhysicalGroup(3, [frag[1][1][2]], 1)
gmsh.model.addPhysicalGroup(3, [frag[1][2][2]], 2)

# set names for outer and inner parts of the battery
gmsh.model.setPhysicalName(3, 1, "outer_volume")
gmsh.model.setPhysicalName(3, 2, "inner_volume")

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate()
gmsh.write("meshes/3d_battery.msh")

gmsh.fltk.run()
gmsh.finalize()