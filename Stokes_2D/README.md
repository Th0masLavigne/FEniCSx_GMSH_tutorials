# Stokes Problem (2D)

The following directory contains the application of a 2D axi-symmetric Stokes Problem. 

The geometry is created using GMSH. Two different pipelines are proposed. It consists in two rectangles of different heights with an inflow imposed on the left side and a no slip condition for the upper sides.

<p align="center">
  <img src="./Stokes2D.jpg">
</p>

The objective is to find the resulting velocity and pressure.

# Implementation

Given the strong compatibility between GMSH and FEniCSx it is recommended to use GMSH. GMSH also has a python API. The mesh, refinement and marking operations can be proceeded in GMSH and imported in the FEniCSx environment.

## Geomerty creation using GMSG
### From a Sketch

#### Libraries and settings

As for any GMSH API code, one must first import the libraries

```python
import gmsh
import numpy
import sys
```

The geometrical and dimension of the problem can then be specified:
```python
L1 = 2.
H1 = 1.
L2 = 2.
H2 = 0.5
# Dimension of the problem,
gdim = 3
```

The gmsh model is initialised with its internal settings specified:
```python
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
```

#### Geometry

To compute the geometry, we first put the points of the rectangles:
```python
# Left square
A = gmsh.model.occ.addPoint(0, 0, lc)
B = gmsh.model.occ.addPoint(0, H1, lc)
C = gmsh.model.occ.addPoint(L1, H1, lc)
D = gmsh.model.occ.addPoint(L1, 0, lc)
# Right square
E = gmsh.model.occ.addPoint((L1+L2), H2, lc)
F = gmsh.model.occ.addPoint((L1+L2), 0, lc)
G = gmsh.model.occ.addPoint(L1, H2, lc)
```

Then, we link the points:
```python
lAB  = gmsh.model.occ.addLine(A, B)
lBC  = gmsh.model.occ.addLine(B, C)
lCG  = gmsh.model.occ.addLine(C, G)
lGD  = gmsh.model.occ.addLine(G, D)
lDA  = gmsh.model.occ.addLine(D, A)
# 
lDF = gmsh.model.occ.addLine(D, F)
lFE = gmsh.model.occ.addLine(F,E)
lEG = gmsh.model.occ.addLine(E,G)
```

The surfaces are defined from a curve loop:
```python
cl1 = gmsh.model.occ.addCurveLoop([lAB, lBC, lCG, lGD, lDA])
s1  = gmsh.model.occ.addPlaneSurface([1], cl1)
gmsh.model.occ.synchronize()
cl2 = gmsh.model.occ.addCurveLoop([lGD, lDF, lFE, lEG])
s2  = gmsh.model.occ.addPlaneSurface([2], cl2)
gmsh.model.occ.synchronize()
```

As a safeguard, we ensure removing all duplicates and export the geometry fot control:
```python
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
# 
gmsh.model.occ.synchronize()
gmsh.write('2D_Stokes.geo_unrolled')
```

#### Marking

GMSH creates the mesh for physical groups. Each of these groups are marked. All lines, surfaces and volumes of the model can be listed with:

```python
lines, surfaces, volumes = [gmsh.model.getEntities(d) for d in [1, 2, 3]]
```

It is then required to create empty lists and tag values:

```python
left, top_left, middle_up, top_right, right, bottom = [], [], [], [], [], []
tag_left, tag_top_left, tag_middle_up, tag_top_right, tag_right, tag_bottom = 1, 2, 3, 4, 5, 6
left_surf, right_surf = [], []
tag_left_surf, tag_right_surf = 10, 20
```

The lists can be automatically filled by identification of the faces and volumes based on their center of mass position:
```python
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
```

Alternatively they can be hand filled using the geo_unrolled filed and visualizing in the GMSH GUI.

To assign the surface tags, we run the following:
```python
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
```

Similarly for the volumes:
```python
for surface in surfaces:
  center_of_mass = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
  if center_of_mass[0]<L1:
    left_surf.append(surface[1])
  else:
    right_surf.append(surface[1])
# 
gmsh.model.addPhysicalGroup(gdim-1, left_surf, tag_left_surf)
gmsh.model.setPhysicalName(gdim-1, tag_left_surf, 'left')
# 
gmsh.model.addPhysicalGroup(gdim-1, right_surf, tag_right_surf)
gmsh.model.setPhysicalName(gdim-1, tag_right_surf, 'right')
```

Once again it is recommended to check the geometry identification:
```python
gmsh.model.occ.synchronize()
gmsh.write('2D_Stokes_marked.geo_unrolled')
```

The mesh is generated and exported:
```python
gmsh.model.mesh.generate(gdim-1)
gmsh.write("2D_Stokes_mesh.msh")
```

Executing the following command at the end run the GMSH Gui for visualization before finalizing the model:
```python
if 'close' not in sys.argv:
    gmsh.fltk.run()
```
```python
gmsh.finalize()
```

### From elementary geometries
When you have elementary entities defining your domain, it is often easier to use them and use boolean operations (cut,fragment,fuse,*etc.*). In the present case the Geometry could be created simply with two rectangles. 

The beginning is the same:

```python
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
```

Then, two rectangles are added and our geometry is finished:
```python
# gmsh.model.occ.addRectangle(x, y, z, dx, dy, tag=-1, roundedRadius=0.)
s1 = gmsh.model.occ.addRectangle(A[0], A[1], 0, L1, H1, tag=-1)
s2 = gmsh.model.occ.addRectangle(D[0], D[1], 0, L2, H2, tag=-1)
```
The duplicates are removed and a check file is generated.
```python
# Remove duplicate entities and synchronize
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
# 
gmsh.model.occ.synchronize()
gmsh.write('2D_Stokes.geo_unrolled')
```

All the following steps (marking, meshing) are the same as previously:

```python
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
```
```python
if 'close' not in sys.argv:
    gmsh.fltk.run()
```
```python
gmsh.finalize()
```

## Finite Element Computation

### Libraries
Computing the Finite element problem within FEniCSx in python requires to load the libraries:

```python
import dolfinx
import ufl
import basix
import mpi4py
import petsc4py
from dolfinx.fem.petsc import LinearProblem
```

One can assess the version of FEniCSx with the following:
```python
print("Dolfinx version is:",dolfinx.__version__)
```

### Mesh Loading
We load the mesh, facet and cell tags from the msh file created:

```python
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh('./2D_Stokes_mesh.msh', mpi4py.MPI.COMM_WORLD, 0, gdim=2)
```

### Function spaces, Functions and expressions

For stability reasons, we create a mixed element P2vxP1 for (u,p) to solve the Stokes problem:
```python
# Finite Element 
P1       = basix.ufl.element("P", mesh.topology.cell_name(), degree=1)
# Vector Element
P1_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=1, shape=(mesh.topology.dim,))
P2_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=2, shape=(mesh.topology.dim,))
# Mixed element
MxE      = basix.ufl.mixed_element([P1,P2_v])
# Tensor Element to compute the stress
tensor_elem  = basix.ufl.element("P", mesh.topology.cell_name(), degree=1, shape=(3,3))
```

The associated function spaces with the element types are:
```python
P1v_space = dolfinx.fem.functionspace(mesh, P1_v)
CHS       = dolfinx.fem.functionspace(mesh, MxE)
tensor_space = dolfinx.fem.functionspace(mesh, tensor_elem)
```

One can then specify the functions and if needed their expression for interpolation:
```python
u_export      = dolfinx.fem.Function(P1v_space)
u_export.name = "u"
# 
sol    = ufl.TrialFunction(CHS)
q, w   = ufl.TestFunctions(CHS)
# Solution vector
p, u   = ufl.split(sol)
# 
# Definition of the normal vector
n      = ufl.FacetNormal(mesh)
```

The operators are also defined at this point:
```python
# Specify the desired quadrature degree
q_deg = 4
# Redefinition dx and ds
dx    = ufl.Measure('dx', metadata={"quadrature_degree":q_deg}, subdomain_data=cell_tag, domain=mesh)
ds    = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
```

### Boundary conditions