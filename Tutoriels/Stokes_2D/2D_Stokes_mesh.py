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
A = gmsh.model.occ.addPoint(0, 0, lc)
B = gmsh.model.occ.addPoint(0, H1, lc)
C = gmsh.model.occ.addPoint(L1, H1, lc)
D = gmsh.model.occ.addPoint(L1, 0, lc)
# Right square
E = gmsh.model.occ.addPoint((L1+L2), H2, lc)
F = gmsh.model.occ.addPoint((L1+L2), 0, lc)
G = gmsh.model.occ.addPoint(L1, H2, lc)
# 
#Define the lines
# 
lAB  = gmsh.model.occ.addLine(A, B)
lBC  = gmsh.model.occ.addLine(B, C)
lCG  = gmsh.model.occ.addLine(C, G)
lGD  = gmsh.model.occ.addLine(G, D)
lDA  = gmsh.model.occ.addLine(D, A)
# 
lDF = gmsh.model.occ.addLine(D, F)
lFE = gmsh.model.occ.addLine(F,E)
lEG = gmsh.model.occ.addLine(E,G)
# 
# Create the curve loops and surfaces
cl1 = gmsh.model.occ.addCurveLoop([lAB, lBC, lCG, lGD, lDA])
s1  = gmsh.model.occ.addPlaneSurface([1], cl1)
gmsh.model.occ.synchronize()
cl2 = gmsh.model.occ.addCurveLoop([lGD, lDF, lFE, lEG])
s2  = gmsh.model.occ.addPlaneSurface([2], cl2)
gmsh.model.occ.synchronize()
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
for surface in surfaces:
	center_of_mass = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
	if center_of_mass[0]<L1:
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
gmsh.write('2D_Stokes_marked.geo_unrolled')
# 
gmsh.model.mesh.generate(gdim)
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