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
L1 = 2.
H1 = 1.
L2 = 2.
H2 = 0.5
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
gmsh.model.add("2D_Stokes")
#Characteristic length
lc = (L1+L2)/60
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
#Add the points
# Left square
A = [0, 0]
C = [L1, H1]
D = [L1, 0]
# Right square
E = [(L1+L2), H2]
# 
# gmsh.model.occ.addRectangle(x, y, z, dx, dy, tag=-1, roundedRadius=0.)
s1 = gmsh.model.occ.addRectangle(A[0], A[1], 0, L1, H1, tag=-1)
s2 = gmsh.model.occ.addRectangle(D[0], D[1], 0, L2, H2, tag=-1)
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
gmsh.write('2D_Stokes.geo_unrolled')
# 
#----------------------------------------------------------------------
# Create physical group for mesh generation and tagging
#----------------------------------------------------------------------
# 
lines, surfaces, volumes = [gmsh.model.getEntities(d) for d in [1, 2, 3]]
# 
left, top_left, middle_up, top_right, right, bottom = [], [], [], [], [], []
tag_left, tag_top_left, tag_middle_up, tag_top_right, tag_right, tag_bottom = 1, 2, 3, 4, 5, 6
left_surf, right_surf = [], []
tag_left_surf, tag_right_surf = 10, 20
# 
for line in lines:
	center_of_mass = gmsh.model.occ.getCenterOfMass(line[0], line[1])
	if numpy.isclose(center_of_mass[0],0):
		left.append(line[1])
	elif numpy.isclose(center_of_mass[1],H1):
		top_left.append(line[1])
	elif numpy.isclose(center_of_mass[1],(H2+H1)/2):
		middle_up.append(line[1])
	elif numpy.isclose(center_of_mass[1],H2):
		top_right.append(line[1])
	elif numpy.isclose(center_of_mass[0],L1+L2):
		right.append(line[1])
	elif numpy.isclose(center_of_mass[1],0):
		bottom.append(line[1])
# 
gmsh.model.addPhysicalGroup(gdim-2, left, tag_left)
gmsh.model.setPhysicalName(gdim-2, tag_left, 'Left')
# 
gmsh.model.addPhysicalGroup(gdim-2, top_left, tag_top_left)
gmsh.model.setPhysicalName(gdim-2, tag_top_left, 'Top_left')
# 
gmsh.model.addPhysicalGroup(gdim-2, middle_up, tag_middle_up)
gmsh.model.setPhysicalName(gdim-2, tag_middle_up, 'Middle_up')
# 
gmsh.model.addPhysicalGroup(gdim-2, top_right, tag_top_right)
gmsh.model.setPhysicalName(gdim-2, tag_top_right, 'Top_right')
# 
gmsh.model.addPhysicalGroup(gdim-2, right, tag_right)
gmsh.model.setPhysicalName(gdim-2, tag_right, 'Right')
# 
gmsh.model.addPhysicalGroup(gdim-2, bottom, tag_bottom)
gmsh.model.setPhysicalName(gdim-2, tag_bottom, 'Bottom')
# 
for surface in surfaces:
	center_of_mass = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
	if center_of_mass[0]<L1:
		left_surf.append(surface[1])
	else:
		right_surf.append(surface[1])
# 
print(surfaces)
gmsh.model.addPhysicalGroup(gdim-1, left_surf, tag_left_surf)
gmsh.model.setPhysicalName(gdim-1, tag_left_surf, 'left')
# 
gmsh.model.addPhysicalGroup(gdim-1, right_surf, tag_right_surf)
gmsh.model.setPhysicalName(gdim-1, tag_right_surf, 'right')
#----------------------------------------------------------------------
# Export the geometry with the tags for control
#----------------------------------------------------------------------
# 
gmsh.model.occ.synchronize()
gmsh.write('2D_Stokes_marked.geo_unrolled')
# 
gmsh.model.mesh.generate(gdim-1)
gmsh.write("2D_Stokes_mesh.msh")
# 
#----------------------------------------------------------------------
# Observe in gmsh
#----------------------------------------------------------------------
# 
if 'close' not in sys.argv:
	gmsh.fltk.run()
# 
gmsh.finalize()
# EoF