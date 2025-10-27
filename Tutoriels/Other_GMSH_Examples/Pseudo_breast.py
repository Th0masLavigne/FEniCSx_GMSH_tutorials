import numpy as np
import math
import os
import sys
import gmsh

# Display the objective
def print_objective():
	print("""
		The objective of this code is to mesh the following revolution geometry (axis of revolution is ':') with tetrahedrons (! No prisms)

		More specifically, this is a 2-part mesh (Region + skin of depth epsilon)

		:
		:
		:B2(6)
		+---   C2 (7)
		|   ---+--
		+---      --- 
		|B(2)---+-   ----
		|       C(3)--   ---
		|             --    --
		|                --  --
		|                   -   --
		|                    -    -
		|                     -    -  D2(8)
		|                      -    +------
		|                       +-------   -------
		|                       D(4)    --------   ------
		|                                       ------  --
		+----------------------------------------------+--+
		:A(1)                                          E(5) E2 (9)
		:
		:
		:


		Then the objective is to extract it as 2 separate submeshes WITH coinciding nodes on the interface

		First Mesh: The inner part

		+---       
		|B(2)---+-   
		|       C(3)--   
		|             --    
		|                --  
		|                   -   
		|                    -    
		|                     -      
		|                      -    
		|                       +-------   
		|                       D(4)    --------   
		|                                       ------  
		+----------------------------------------------+
		A(1)                                          E(5) 



		Second Mesh: the skin

		:B2(6)
		+---   C2 (7)
		|   ---+--
		+---      --- 
		B(2)---+-   ----
		     C(3)--   ---
		           --    --
		             --  --
		               -   --
		                -    -
		                 -    -  D2(8)
		                  -    +------
		                   +-------   -------
		                           --------   ------
		                                    ------  --
		                                           +--+
		                                           E(5) E2 (9)
		
		""")
	return None

	

# Function returning the center of a circle passing by 3 points P1, P2, P3
def find_center(P1,P2,P3):
	a = P3[1]**2-P2[1]**2+P3[2]**2-P2[2]**2
	b = P3[2]-P2[2]
	c = P2[1]**2-P1[1]**2+P2[2]**2-P1[2]**2
	d = P2[2]-P1[2]
	e = P2[1]-P1[1]
	f = P3[1]-P2[1]
	g = P3[2]-P2[2]
	x = (a/(2*b)-c/(2*d))/(e/d-f/g);
	y = -e/d*x+c/(2*d);
	return x, y

# Function returning points coordinates on a circle (modifies y of target point)
def find_y_point_of_circle(A,O,B):
	R=np.sqrt((A[1]-O[1])**2+(A[2]-O[2])**2)
	yB = O[1]+np.sqrt(R**2-(B[2]-O[2])**2)
	B[1] = yB
	return B

# Function returning points coordinates on a circle (modifies z of target point)
def find_z_point_of_circle(A,O,B):
	R=np.sqrt((A[1]-O[1])**2+(A[2]-O[2])**2)
	zB = O[2]+np.sqrt(R**2-(B[1]-O[1])**2)
	B[2] = zB
	return B

# Function to create the geometry and mesh
def create_geometry_and_mesh(A,B,B2,C,C2,D,D2,E,E2,O1,O2,O3,Mesh_name,output):
	gmsh.initialize()
	#
	#Characteristic length
	lc = 0.005
	gmsh.model.occ.synchronize()
	gmsh.option.setNumber("General.Terminal",1)
	gmsh.option.setNumber("Mesh.Optimize", True)
	gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
	#gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
	#gmsh.option.setNumber("Mesh.MeshSizeMax", 40*lc)
	gmsh.model.occ.synchronize()
	# gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
	# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0.002)
	# gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
	# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5)
	# Define points coordinates for the geometry
	#   B2 (6)
	#	+		C2 (7)
	#	+   	+
	#	B (2)  +
	#			C (3)
	#	
	#	
	#	
	#	
	#							  D2 (8)
	#							  +
	#							+
	#							D (4)
	#	
	#	+                                    + +
	#	A (1)									 E(5) E2 (9)
	#	
	#
	# Circle arc : B,C and B2,C2 => center O1
	# Circle arc : C,D and C2,D2 => center O2
	# Circle arc : D,E and D2,E2 => center O3
	# Skin depth is epsilon
	# Clear all models and create a new one
	gmsh.clear()
	gmsh.model.add("Skinned_Breast_Phantom")
	#gmsh.model.occ.healShapes(dimTags = [], tolerance = 1e-5, fixDegenerated = True, fixSmallEdges = True, fixSmallFaces = True, sewFaces = True, makeSolids = True )
	#
	#Add the points
	Ag  = gmsh.model.occ.addPoint(A[0],  A[1],  A[2],  400*lc)
	Bg  = gmsh.model.occ.addPoint(B[0],  B[1],  B[2],  100*lc)
	Cg  = gmsh.model.occ.addPoint(C[0],  C[1],  C[2],  lc)
	Dg  = gmsh.model.occ.addPoint(D[0],  D[1],  D[2],  lc)
	Eg  = gmsh.model.occ.addPoint(E[0],  E[1],  E[2],  100*lc)
	B2g = gmsh.model.occ.addPoint(B2[0], B2[1], B2[2], lc)
	C2g = gmsh.model.occ.addPoint(C2[0], C2[1], C2[2], lc)
	D2g = gmsh.model.occ.addPoint(D2[0], D2[1], D2[2], lc)
	E2g = gmsh.model.occ.addPoint(E2[0], E2[1], E2[2], lc)
	#
	O1g = gmsh.model.occ.addPoint(O1[0], O1[1], O1[2], lc)
	O2g = gmsh.model.occ.addPoint(O2[0], O2[1], O2[2], lc)
	O3g = gmsh.model.occ.addPoint(O3[0], O3[1], O3[2], lc)
	#
	#Define the lines
	lAB  = gmsh.model.occ.addLine(Ag, Bg)
	lAE  = gmsh.model.occ.addLine(Ag, Eg)
	lEE2 = gmsh.model.occ.addLine(Eg, E2g)
	lBB2 = gmsh.model.occ.addLine(Bg, B2g)
	#Define the circle arcs (start point, center, end point)
	lBC    = gmsh.model.occ.addCircleArc(Bg,O1g,Cg)
	lB2C2  = gmsh.model.occ.addCircleArc(B2g,O1g,C2g)
	lCD    = gmsh.model.occ.addCircleArc(Cg,O2g,Dg)
	lC2D2  = gmsh.model.occ.addCircleArc(C2g,O2g,D2g)
	lDE    = gmsh.model.occ.addCircleArc(Dg,O3g,Eg)
	lD2E2  = gmsh.model.occ.addCircleArc(D2g,O3g,E2g)
	#
	# Create the curve loops and surfaces
	cl1 = gmsh.model.occ.addCurveLoop([lAB, lBC, lCD, lDE, -lAE])
	s1  = gmsh.model.occ.addPlaneSurface([1], cl1)
	gmsh.model.occ.synchronize()
	#
	#
	cl2 = gmsh.model.occ.addCurveLoop([lBB2, lB2C2, lC2D2, lD2E2, -lEE2, -lDE, -lCD, -lBC])
	s2 = gmsh.model.occ.addPlaneSurface([2], cl2)
	gmsh.model.occ.synchronize()
	#
	import os
	gmsh.write('Geom_reelle_8EP.geo_unrolled')
	#
	# Make the revolution geometry
	ov = gmsh.model.occ.revolve([(2,s1)], 0, 0, 0, 0, 0, 1, 2*math.pi,[1])
	ov2 = gmsh.model.occ.revolve([(2,s2)], 0, 0, 0, 0, 0, 1, 2*math.pi,[1])
	# Synchronize
	gmsh.model.occ.synchronize()
	#
	# Remove duplicate entities and synchronize
	gmsh.model.occ.removeAllDuplicates()
	gmsh.model.occ.synchronize()
	surfaces, volumes = [gmsh.model.getEntities(d) for d in [ 2, 3]]
	print(surfaces)
	print(volumes)
	if output == 0:
		tdim =3 
		gmsh.model.addPhysicalGroup(tdim, [1], 1)
		gmsh.model.setPhysicalName(tdim, 1, 'Tissue')
	elif output == 1:
		tdim =3
		gmsh.model.addPhysicalGroup(tdim, [2], 2)
		gmsh.model.setPhysicalName(tdim, 2, 'Skin')
	elif output == 2:
		tdim =3 
		gmsh.model.addPhysicalGroup(tdim, [1], 1)
		gmsh.model.setPhysicalName(tdim, 1, 'Tissue')
		gmsh.model.addPhysicalGroup(tdim, [2], 2)
		gmsh.model.setPhysicalName(tdim, 2, 'Skin')
		tdim =2 
		gmsh.model.addPhysicalGroup(tdim, [3,4,5], 5)
		gmsh.model.setPhysicalName(tdim, 5, 'Interface')
		gmsh.model.addPhysicalGroup(tdim, [7,8,9], 6)
		gmsh.model.setPhysicalName(tdim, 6, 'Skin Surface')
	elif output == 2:
		tdim =3 
		gmsh.model.addPhysicalGroup(tdim, [1], 1)
		gmsh.model.setPhysicalName(tdim, 1, 'Tissue')
		gmsh.model.addPhysicalGroup(tdim, [2], 2)
		gmsh.model.setPhysicalName(tdim, 2, 'Skin')
	#
	# We finally generate and save the mesh:
	gmsh.model.mesh.generate(3)
	gmsh.write(Mesh_name+".msh")
	gmsh.finalize()
	return gmsh.model


if __name__ == "__main__":

	print_objective()

	epsilon = 0.005

	A  = [0,	0,							0]
	B  = [0,	0, 							0.08]
	B2 = [0, 	0,							0.08+epsilon]
	C  = [0,	0.0075, 					0.075]
	C2 = [0, 	0.0075+ epsilon/np.sqrt(2), 0.075+epsilon/np.sqrt(2)]
	D  = [0,	0.08,						0.01]
	D2 = [0, 	0.08+epsilon/np.sqrt(2), 	0.01+epsilon/np.sqrt(2)]
	E  = [0,	0.14,						0]
	E2 = [0, 	0.14+epsilon, 				0]

	O1 = [0,	-0.0002,				0.0721]
	O2 = [0,	0.0039,				-0.0014]
	O3 = [0,	0.0049,					-0.6385]

	print("Beware to update the positions in the order: E/E2, C/C2, B/B2 \n")

	E  = find_y_point_of_circle(D,O3,E)
	E2 = find_y_point_of_circle(D2,O3,E2)
	C  = find_y_point_of_circle(D,O2,C)
	C2 = find_y_point_of_circle(D2,O2,C2)
	B  = find_z_point_of_circle(C,O1,B)
	B2 = find_z_point_of_circle(C2,O1,B2)


	# 2 = all
	part = 2;
	Mesh_name = "All"
	gmsh_model=create_geometry_and_mesh(A,B,B2,C,C2,D,D2,E,E2,O1,O2,O3,Mesh_name,part)

	# 3 = all volume tetra
	part = 3;
	Mesh_name = "All_vol"
	gmsh_model=create_geometry_and_mesh(A,B,B2,C,C2,D,D2,E,E2,O1,O2,O3,Mesh_name,part)

	# 0 = tissue
	part = 0;
	Mesh_name = "Tissue"
	gmsh_model=create_geometry_and_mesh(A,B,B2,C,C2,D,D2,E,E2,O1,O2,O3,Mesh_name,part)

	# 1 = skin
	part = 1;
	Mesh_name = "Skin"
	gmsh_model=create_geometry_and_mesh(A,B,B2,C,C2,D,D2,E,E2,O1,O2,O3,Mesh_name,part)
