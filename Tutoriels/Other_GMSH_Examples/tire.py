import numpy as np
import math
import os
import sys
import gmsh

# Display the objective
def print_objective():
	print("""
		The objective of this code is to mesh the following revolution geometry (axis of revolution is ':') 





                                           G(7)
									   ----+----
								    ---          ---
		I(9)                 H(8)---      O(11)    --  F(6)
		+------------------------+         +          --+------------------------+E(5)
		|                        |                      |                        |
		+------------------------+----------------------+------------------------+
		A(1)                     B(2)                   C(3)                      D(4)

		
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
def create_geometry_and_mesh(A,B,C,D,E,F,H,I,O,Mesh_name,output):
	gmsh.initialize()
	#
	#Characteristic length
	#lc = 0.02
	lc = 0.05
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
  	#                                  G(7)
	# 							   ----+----
	# 						    ---          ---
	# I(9)                 H(8)---      O(11)    --  F(6)
	# +------------------------+         +          --+------------------------+E(5)
	# |                        |                      |                        |
	# +------------------------+----------------------+------------------------+
	# A(1)                     B(2)                   C(3)                      D(4)


	#	
	# Clear all models and create a new one
	gmsh.clear()
	gmsh.model.add("Contact_Geometry")
	#gmsh.model.occ.healShapes(dimTags = [], tolerance = 1e-5, fixDegenerated = True, fixSmallEdges = True, fixSmallFaces = True, sewFaces = True, makeSolids = True )
	#
	#Add the points
	Ag  = gmsh.model.occ.addPoint(A[0],  A[1],  A[2],  lc)
	Bg  = gmsh.model.occ.addPoint(B[0],  B[1],  B[2],  lc)
	Cg  = gmsh.model.occ.addPoint(C[0],  C[1],  C[2],  lc)
	Dg  = gmsh.model.occ.addPoint(D[0],  D[1],  D[2],  lc)
	Eg  = gmsh.model.occ.addPoint(E[0],  E[1],  E[2],  lc)
	Fg  = gmsh.model.occ.addPoint(F[0],  F[1],  F[2],  lc)
	Gg  = gmsh.model.occ.addPoint(G[0],  G[1],  G[2],  lc)
	Hg  = gmsh.model.occ.addPoint(H[0],  H[1],  H[2],  lc)
	Ig  = gmsh.model.occ.addPoint(I[0],  I[1],  I[2],  lc)
	Og  = gmsh.model.occ.addPoint(O[0],  O[1],  O[2],  lc)
	#
	#Define the lines
	lAB  = gmsh.model.occ.addLine(Ag, Bg)
	lBH  = gmsh.model.occ.addLine(Bg, Hg)
	lHI  = gmsh.model.occ.addLine(Hg, Ig)
	lAI  = gmsh.model.occ.addLine(Ag, Ig)
	lFC  = gmsh.model.occ.addLine(Fg, Cg)
	lCD  = gmsh.model.occ.addLine(Cg, Dg)
	lDE  = gmsh.model.occ.addLine(Dg, Eg)
	lEF  = gmsh.model.occ.addLine(Eg, Fg)
	lBC = gmsh.model.occ.addLine(Bg, Cg)
	#Define the circle arcs (start point, center, end point)
	lFG    = gmsh.model.occ.addCircleArc(Fg,Og,Gg)
	lGH    = gmsh.model.occ.addCircleArc(Gg,Og,Hg)
	#
	# Create the curve loops and surfaces
	cl1 = gmsh.model.occ.addCurveLoop([lAB, lBH, lHI,-lAI])
	s1  = gmsh.model.occ.addPlaneSurface([1], cl1)
	gmsh.model.occ.synchronize()
	#
	#
	#
	cl2 = gmsh.model.occ.addCurveLoop([lCD, lDE, lEF, lFC])
	s2 = gmsh.model.occ.addPlaneSurface([2], cl2)
	gmsh.model.occ.synchronize()
	#
	cl3 = gmsh.model.occ.addCurveLoop([lBC, -lFC, lFG, lGH, -lBH])
	s3 = gmsh.model.occ.addPlaneSurface([3], cl3)
	gmsh.model.occ.synchronize()
	#
	import os
	gmsh.write('Contact_Geo_8EP.geo_unrolled')
	#
	# Make the revolution geometry
	ov = gmsh.model.occ.revolve([(2,s1)], 0, 0, 0, 0, 0, 1, 2*math.pi,[1])
	ov2 = gmsh.model.occ.revolve([(2,s2)], 0, 0, 0, 0, 0, 1, 2*math.pi,[1])
	ov3 = gmsh.model.occ.revolve([(2,s3)], 0, 0, 0, 0, 0, 1, 2*math.pi,[1])
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
		gmsh.model.addPhysicalGroup(tdim, [1,2], 1)
		gmsh.model.setPhysicalName(tdim, 1, 'Cylinder')
	elif output == 1:
		tdim =3
		gmsh.model.addPhysicalGroup(tdim, [3], 2)
		gmsh.model.setPhysicalName(tdim, 2, 'Sphere')
	#
	# We finally generate and save the mesh:
	gmsh.model.mesh.generate(3)
	gmsh.write(Mesh_name+".msh")
	gmsh.finalize()
	return gmsh.model


if __name__ == "__main__":

	print_objective()
	A  = [0, 0.0, 0.0]
	B  = [0, 0.0, 0.05]
	C  = [0, 0.0, 0.25]
	D  = [0, 0.0, 0.3]
	E  = [0, 0.05, 0.3]
	F  = [0, 0.05, 0.25]
	G  = [0, 0.15, 0.15]
	H  = [0, 0.05, 0.05]
	I  = [0, 0.05, 0.0]
	O  = [0, 0.05, 0.15]

	# print(Bpos)
	# Bpos  = find_z_point_of_circle(Cpos,Epos,Bpos)
	# print(Bpos)

	Mesh_name = "Sphere"
	gmsh_model=create_geometry_and_mesh(A,B,C,D,E,F,H,I,O,Mesh_name,1)
	Mesh_name = "Cyl"
	gmsh_model=create_geometry_and_mesh(A,B,C,D,E,F,H,I,O,Mesh_name,0)

