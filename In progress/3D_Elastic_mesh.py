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
# Dimension of the problem,
gdim = 3
#----------------------------------------------------------------------
# 
#----------------------------------------------------------------------
# Set options
#----------------------------------------------------------------------
# 
gmsh.initialize()
gmsh.clear()
gmsh.model.add("3D_Elastic")
#Characteristic length
lc = 40/20
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
# 
s1 = gmsh.model.occ.addBox(0, 0, 0, 20, 1, 1, tag=-1)
s2 = gmsh.model.occ.addBox(20, 0, 0, 20, 1, 1, tag=-1)
# s1 = gmsh.model.occ.addRectangle(0, 0, 0, 20, 5, tag=-1)
# s2 = gmsh.model.occ.addRectangle(20, 0, 0, 20, 5, tag=-1)
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
gmsh.write('3D_elastic.geo_unrolled')
# 
#----------------------------------------------------------------------
# Create physical group for mesh generation and tagging
#----------------------------------------------------------------------
# 
lines, surfaces, volumes = [gmsh.model.getEntities(d) for d in [1, 2, 3]]
# 
left, top, right, bottom = [], [], [], []
tag_left, tag_top, tag_right, tag_bottom = 1, 2, 3, 4
left_surf, right_surf = [], []
tag_left_surf, tag_right_surf = 10, 20
# 
for border in surfaces:
	center_of_mass = gmsh.model.occ.getCenterOfMass(border[0], border[1])
	if numpy.isclose(center_of_mass[0],0):
		left.append(border[1])
	elif numpy.isclose(center_of_mass[2],1):
		top.append(border[1])
	elif numpy.isclose(center_of_mass[0],40):
		right.append(border[1])
	elif numpy.isclose(center_of_mass[2],0):
		bottom.append(border[1])
# 
gmsh.model.addPhysicalGroup(gdim-1, left, tag_left)
gmsh.model.setPhysicalName(gdim-1, tag_left, 'Left')
# 
gmsh.model.addPhysicalGroup(gdim-1, top, tag_top)
gmsh.model.setPhysicalName(gdim-1, tag_top, 'Top')
# 
gmsh.model.addPhysicalGroup(gdim-1, right, tag_right)
gmsh.model.setPhysicalName(gdim-1, tag_right, 'Right')
# 
gmsh.model.addPhysicalGroup(gdim-1, bottom, tag_bottom)
gmsh.model.setPhysicalName(gdim-1, tag_bottom, 'Bottom')
# 
for surface in surfaces:
	center_of_mass = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
	if center_of_mass[0]<20:
		left_surf.append(surface[1])
	else:
		right_surf.append(surface[1])
# 
gmsh.model.addPhysicalGroup(gdim, left_surf, tag_left_surf)
gmsh.model.setPhysicalName(gdim, tag_left_surf, 'left')
# 
gmsh.model.addPhysicalGroup(gdim, right_surf, tag_right_surf)
gmsh.model.setPhysicalName(gdim, tag_right_surf, 'right')
#----------------------------------------------------------------------
# Export the geometry with the tags for control
#----------------------------------------------------------------------
# 
gmsh.model.occ.synchronize()
gmsh.write('3D_elastic.geo_unrolled')
# 
gmsh.model.mesh.generate(gdim)
gmsh.write("3D_elastic_mesh.msh")
# 
#----------------------------------------------------------------------
# Observe in gmsh
#----------------------------------------------------------------------
# 
if 'close' not in sys.argv:
	gmsh.fltk.run()
# 
gmsh.finalize()