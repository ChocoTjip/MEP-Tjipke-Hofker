using FerriteGmsh


gmsh.initialize()

gmsh.model.add("battery_2d_corner")

h1 = .1
w1 = .3
i1 = .04

lc1 = 1e-2
lc2 = 2e-2


#redefine namespace to make it easier
factory = gmsh.model.occ

# define corners
p1 = factory.addPoint(-w1/2, -h1/2, 0, lc1)
p2 = factory.addPoint(-w1/2,  h1/2, 0, lc1)
p3 = factory.addPoint( w1/2,  h1/2, 0, lc1)
p4 = factory.addPoint( w1/2, -h1/2, 0, lc1)

p5 = factory.addPoint(-w2/2, -h2/2, 0, lc2)
p6 = factory.addPoint(-w2/2,  h2/2, 0, lc2)
p7 = factory.addPoint( w2/2,  h2/2, 0, lc2)
p8 = factory.addPoint( w2/2, -h2/2, 0, lc2)

p9  = factory.addPoint(-w1/2, -i1/2, 0, lc3)
p10 = factory.addPoint(-w2/2, -i1/2, 0, lc3)
p11 = factory.addPoint(-w2/2,  i1/2, 0, lc3)
p12 = factory.addPoint(-w1/2,  i1/2, 0, lc3)

# l1 = factory.addLine(p1, p2)
l2 = factory.addLine(p2, p3)
l3 = factory.addLine(p3, p4)
l4 = factory.addLine(p4, p1)

l5 = factory.addLine(p5, p6)
l6 = factory.addLine(p6, p7)
l7 = factory.addLine(p7, p8)
l8 = factory.addLine(p8, p5)

l9  = factory.addLine(p1 , p9)
l10 = factory.addLine(p9 , p10)
l11 = factory.addLine(p10, p11)
l12 = factory.addLine(p11, p12)
l13 = factory.addLine(p12, p2)

loop1 = factory.addCurveLoop([l9, l10, l11, l12, l13, l2, l3, l4])
loop2 = factory.addCurveLoop([l5, l6, l7, l8])

inner = factory.addPlaneSurface([loop1])
outer = factory.addPlaneSurface([loop2])

frag = factory.fragment([(2, outer)], [(2, inner)])

println(gmsh.model.getEntities(2))

# frag[1] gives new entities from the fragment function
# frag[1][1] gives first (dim, tag)-tuple of the old fragment function
# frag[1][1][2] gives tag of outer surface
gmsh.model.addPhysicalGroup(2, [frag[1][1][2]], 1)
gmsh.model.addPhysicalGroup(2, [frag[1][2][2]], 2)

# set names for outer and inner parts of the battery
gmsh.model.setPhysicalName(2, 1, "outer_surface")
gmsh.model.setPhysicalName(2, 2, "inner_surface")

factory.synchronize()
gmsh.model.mesh.generate()
gmsh.write("meshes/2d_battery.msh")
# gmsh.option.setNumber("Gui.PhysicalGroups", 1)

gmsh.fltk.run()

gmsh.finalize()