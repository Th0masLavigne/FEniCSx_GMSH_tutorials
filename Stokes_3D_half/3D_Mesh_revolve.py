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
gmsh.model.add("3D_Stokes")
#Characteristic length
lc = (L1+L2)/30
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
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0.5)
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
gmsh.model.occ.synchronize()
# 
gmsh.model.occ.revolve([(2,s1),(2,s2)], 0,0,0, 1, 0, 0, -numpy.pi/2)
gmsh.model.occ.synchronize()
# Remove duplicate entities and synchronize
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
# 
#----------------------------------------------------------------------
# Export the geometry for control
#----------------------------------------------------------------------
# 
gmsh.model.occ.synchronize()
gmsh.write('3D_Stokes.geo_unrolled')
# 
#----------------------------------------------------------------------
# Create physical group for mesh generation and tagging
#----------------------------------------------------------------------
# 
lines, surfaces, volumes = [gmsh.model.getEntities(d) for d in [1, 2, 3]]
boundaries = gmsh.model.getBoundary(volumes, oriented=False)
# 
left, top_left, middle_up, top_right, right, bottom, front = [], [], [], [], [], [], []
tag_left, tag_top_left, tag_middle_up, tag_top_right, tag_right, tag_bottom, tag_front = 1, 2, 3, 4, 5, 6,7
left_vol, right_vol = [], []
tag_left_vol, tag_right_vol = 10, 20
# 
for boundary in boundaries:
	center_of_mass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
	if numpy.isclose(center_of_mass[0],0):
		left.append(boundary[1])
	elif (center_of_mass[0]>H2) and (center_of_mass[1]>H2):
		top_left.append(boundary[1])
	elif numpy.isclose(center_of_mass[0],L1):
		middle_up.append(boundary[1])
	elif numpy.isclose(center_of_mass[0],L1+L2):
		right.append(boundary[1])
	elif numpy.isclose(center_of_mass[1],0):
		bottom.append(boundary[1])
	elif numpy.isclose(center_of_mass[2],0):
		front.append(boundary[1])
	else:
		top_right.append(boundary[1])
# 
gmsh.model.addPhysicalGroup(gdim-1, left, tag_left)
gmsh.model.setPhysicalName(gdim-1, tag_left, 'Left')
# 
gmsh.model.addPhysicalGroup(gdim-1, top_left, tag_top_left)
gmsh.model.setPhysicalName(gdim-1, tag_top_left, 'Top_left')
# 
gmsh.model.addPhysicalGroup(gdim-1, middle_up, tag_middle_up)
gmsh.model.setPhysicalName(gdim-1, tag_middle_up, 'Middle_up')
# 
gmsh.model.addPhysicalGroup(gdim-1, top_right, tag_top_right)
gmsh.model.setPhysicalName(gdim-1, tag_top_right, 'Top_right')
# 
gmsh.model.addPhysicalGroup(gdim-1, right, tag_right)
gmsh.model.setPhysicalName(gdim-1, tag_right, 'Right')
# 
gmsh.model.addPhysicalGroup(gdim-1, bottom, tag_bottom)
gmsh.model.setPhysicalName(gdim-1, tag_bottom, 'Bottom')
# 
gmsh.model.addPhysicalGroup(gdim-1, front, tag_front)
gmsh.model.setPhysicalName(gdim-1, tag_front, 'Front')
# 
for volume in volumes:
	center_of_mass = gmsh.model.occ.getCenterOfMass(volume[0], volume[1])
	if center_of_mass[0]<L1:
		left_vol.append(volume[1])
	else:
		right_vol.append(volume[1])
# 
gmsh.model.addPhysicalGroup(gdim, left_vol, tag_left_vol)
gmsh.model.setPhysicalName(gdim, tag_left_vol, 'left')
# 
gmsh.model.addPhysicalGroup(gdim, right_vol, tag_right_vol)
gmsh.model.setPhysicalName(gdim, tag_right_vol, 'right')
# 
#----------------------------------------------------------------------
# Export the geometry with the tags for control
#----------------------------------------------------------------------
# 
gmsh.model.occ.synchronize()
gmsh.write('3D_Stokes_marked.geo_unrolled')
# 
# 
gmsh.model.mesh.generate(gdim)
gmsh.write("3D_Stokes_mesh.msh")
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