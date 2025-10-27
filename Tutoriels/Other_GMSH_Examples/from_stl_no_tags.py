# Thomas Lavigne
# 03/08/2024
# 
import math
import gmsh
import sys
# 
# OpenCascade Kernel does not support stls. This is a limitation for boolean operations

# There a few libraries out there that support boolean operations for meshes, you might want to give them a try. Here's a personal list I have of mostly C++ and Python packages from a comment of @gnikit:
#  https://gitlab.onelab.info/gmsh/gmsh/-/issues/2493
# - CGAL: https://github.com/CGAL/cgal
# - VTK-based improved booleans: https://github.com/zippy84/vtkbool
# - MeshLib: https://github.com/MeshInspector/MeshLib
# - mcut: https://github.com/cutdigital/mcut
# - trimesh: https://github.com/mikedh/trimesh
# - vedo: https://github.com/marcomusy/vedo/blob/master/examples/basic/boolean.py
# - PyMeshLab: https://github.com/cnr-isti-vclab/PyMeshLab (and its full application MeshLab: https://www.meshlab.net/)
# - PyMesh: https://github.com/PyMesh/PyMesh
# - Pycork: https://pypi.org/project/pycork/ (and it's C base library cork: https://github.com/gilbo/cork)
# - Blender: https://www.blender.org/ (also has a python module)

gmsh.initialize()
gmsh.clear()
gmsh.option.setNumber("General.Terminal",1)
gmsh.option.setNumber("Mesh.Optimize", True)
gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
# gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01*lc)
# gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
# gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
# print("Min mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMin"))
# print("Max mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMax"))

angle = 10
forceParametrizablePatches = True
includeBoundary = True
curveAngle = 30
# 
gmsh.merge("Sphere.stl")
# 
gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                 forceParametrizablePatches)
# 
surfaces, volumes = [gmsh.model.getEntities(d) for d in [ 2, 3]]
gmsh.model.geo.synchronize()

gmsh.model.mesh.createGeometry()
l = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfaces])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()
tdim =3 
gmsh.model.addPhysicalGroup(tdim, [1], 1)
gmsh.model.setPhysicalName(tdim, 1, 'test')
# 
gmsh.model.mesh.generate(3)
gmsh.write("Mesh_Sphere.msh")
# 
if 'close' not in sys.argv:
    gmsh.fltk.run()
#EoF