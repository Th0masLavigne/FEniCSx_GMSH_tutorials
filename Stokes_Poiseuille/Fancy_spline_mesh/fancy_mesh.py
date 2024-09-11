# Essai spline

# Environnement
import gmsh
import numpy
import sys

L = 10

# def sine superior
def sin_up(x):
	return 5+numpy.sin(numpy.pi*x/L)

# def cosine inferior
def cos_down(x):
	return -5+numpy.cos(numpy.pi*x/L)

x = numpy.linspace(0,10*L,num=50)
x2 = numpy.linspace(0,10*L,num=150)

y_top = sin_up(x)

y_bottom = cos_down(x2)

# Dimension of the problem,
gdim = 2

gmsh.initialize()
gmsh.clear()
gmsh.model.add("test_spline")
#Characteristic length
lc = L/15
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

# For the bottom points
for index in range(len(y_bottom)):
	gmsh.model.occ.addPoint(x2[index], y_bottom[index],0, meshSize=lc, tag=-1)


# For the top points
for index in range(len(y_top)):
	gmsh.model.occ.addPoint(x[index], y_top[index],0, meshSize=lc, tag=-1)


gmsh.model.occ.synchronize()

points_tuples= gmsh.model.getEntities(0)

points_tuples_bottom = points_tuples[:len(y_bottom)]

points_tuples_top = points_tuples[len(y_bottom):]


points_tags_top = [tuples[1] for tuples in points_tuples_top]

points_tags_bottom = []
for tuples in points_tuples_bottom:
	points_tags_bottom.append(tuples[1])


spline_top = gmsh.model.occ.addSpline(points_tags_top)
spline_bottom = gmsh.model.occ.addSpline(points_tags_bottom, tag=600)

# Sides
line_left = gmsh.model.occ.addLine(points_tags_top[0],points_tags_bottom[0])
line_right = gmsh.model.occ.addLine(points_tags_top[-1],points_tags_bottom[-1])

# Closed Contour
cl1 = gmsh.model.occ.addCurveLoop([spline_top, line_right, -spline_bottom, -line_left])

# create surface
s1  = gmsh.model.occ.addPlaneSurface([cl1],tag=1000)

gmsh.model.occ.synchronize()
# 
gmsh.write('test_spline.geo_unrolled')
# 
# 
gmsh.model.addPhysicalGroup(gdim-1, [line_left], 1)
gmsh.model.setPhysicalName(gdim-1, 1, 'Left')
# 
gmsh.model.addPhysicalGroup(gdim-1, [spline_top], 2)
gmsh.model.setPhysicalName(gdim-1, 2, 'Top')
# 
# 
gmsh.model.addPhysicalGroup(gdim-1, [line_right], 3)
gmsh.model.setPhysicalName(gdim-1, 3, 'Right')
# 
gmsh.model.addPhysicalGroup(gdim-1, [spline_bottom], 4)
gmsh.model.setPhysicalName(gdim-1, 4, 'Bottom')

gmsh.model.addPhysicalGroup(gdim, [s1], 10)
gmsh.model.setPhysicalName(gdim, 10, 'surface')

gmsh.model.occ.synchronize()
# 
gmsh.write('test_spline.geo_unrolled')
# 
# 

gmsh.model.mesh.generate(gdim)
gmsh.write("2D_Stokes_fancy_mesh.msh")
gmsh.finalize()