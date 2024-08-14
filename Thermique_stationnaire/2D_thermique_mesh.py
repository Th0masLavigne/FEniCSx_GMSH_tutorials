# Thomas Lavigne
# 03/08/2024
# 
#----------------------------------------------------------------------
# Libraries
#----------------------------------------------------------------------
# 
import gmsh
import numpy
import sys
# 
#----------------------------------------------------------------------
# Geometrical parameters
#----------------------------------------------------------------------
r_c = 0.02
l_d = 0.24
h_d = 0.1
#----------------------------------------------------------------------
# Dimension of the problem,
gdim = 2
#----------------------------------------------------------------------
# 
#----------------------------------------------------------------------
# Set options
#----------------------------------------------------------------------
# 
gmsh.initialize()
gmsh.clear()
gmsh.model.add("2D_Stokes")
#Characteristic length
lc = l_d/60
gmsh.model.occ.synchronize()
gmsh.option.setNumber("General.Terminal",1)
gmsh.option.setNumber("Mesh.Optimize", True)
gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1*lc)
gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
gmsh.model.occ.synchronize()
# gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0.002)
# gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5)
# 
#----------------------------------------------------------------------
# Compute the geometry
#----------------------------------------------------------------------
# 
# gmsh.model.occ.addRectangle(x, y, z, dx, dy, tag=-1, roundedRadius=0.)
s1 = gmsh.model.occ.addRectangle(-0.5*l_d, 0, 0, l_d, h_d, tag=-1)
gmsh.model.occ.synchronize()
s2 = gmsh.model.occ.addDisk(0, 0, 0, r_c,r_c, tag=-1, zAxis=[], xAxis=[])
# gmsh.model.occ.removeAllDuplicates()
# gmsh.model.occ.synchronize()
# 
# 
domain = gmsh.model.occ.cut([(gdim,s1)], [(gdim,s2)], tag=-1, removeObject=True, removeTool=True)
# 
# Remove duplicate entities and synchronize
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
# 
#----------------------------------------------------------------------
# Export the geometry for control
#----------------------------------------------------------------------
# 
gmsh.model.occ.synchronize()
gmsh.write('2D_thermique.geo_unrolled')
# 
#----------------------------------------------------------------------
# Create physical group for mesh generation and tagging
#----------------------------------------------------------------------
# 
lines, surfaces, volumes = [gmsh.model.getEntities(d) for d in [1, 2, 3]]
boundaries = gmsh.model.getBoundary(surfaces, oriented=False)
# 
left, top, right, bottom, circle = [], [], [], [], []
tag_left, tag_top, tag_right, tag_bottom, tag_circle = 1, 2, 3, 4, 5
surf, tag_surf = [], 10
# 
for line in boundaries:
	center_of_mass = gmsh.model.occ.getCenterOfMass(line[0], line[1])
	if numpy.isclose(center_of_mass[0],-0.5*l_d):
		left.append(line[1])
	elif numpy.isclose(center_of_mass[1],h_d):
		top.append(line[1])
	elif numpy.isclose(center_of_mass[0],0.5*l_d):
		right.append(line[1])
	elif numpy.isclose(center_of_mass[1],0):
		bottom.append(line[1])
	else:
		circle.append(line[1])
# 
gmsh.model.addPhysicalGroup(gdim-1, left, tag_left)
gmsh.model.setPhysicalName(gdim-1, tag_left, 'Left')
# 
gmsh.model.addPhysicalGroup(gdim-1, top, tag_top)
gmsh.model.setPhysicalName(gdim-1, tag_top, 'Top')
# 
gmsh.model.addPhysicalGroup(gdim-1, circle, tag_circle)
gmsh.model.setPhysicalName(gdim-1, tag_circle, 'Circle')
# 
gmsh.model.addPhysicalGroup(gdim-1, right, tag_right)
gmsh.model.setPhysicalName(gdim-1, tag_right, 'right')
# 
gmsh.model.addPhysicalGroup(gdim-1, bottom, tag_bottom)
gmsh.model.setPhysicalName(gdim-1, tag_bottom, 'Bottom')
# 
for surface in surfaces:
	center_of_mass = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
	surf.append(surface[1])
# 
gmsh.model.addPhysicalGroup(gdim, surf, tag_surf)
gmsh.model.setPhysicalName(gdim, tag_surf, 'surface')
# 
#----------------------------------------------------------------------
# Export the geometry with the tags for control
#----------------------------------------------------------------------
# 
gmsh.model.occ.synchronize()
gmsh.write('2D_thermique.geo_unrolled')
# 
gmsh.model.mesh.generate(gdim)
gmsh.write("2D_thermique.msh")
# 
#----------------------------------------------------------------------
# Observe in gmsh
#----------------------------------------------------------------------
# 
if 'close' not in sys.argv:
	gmsh.fltk.run()
# 
gmsh.finalize()