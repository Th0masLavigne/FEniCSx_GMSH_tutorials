# Elastic Multimaterial Beam

## Description of the problem

This example aims in providing an example of a 3D beam of dimensions $`40\times5\times5`$ composed of two materials.

The beam is subdivided into two subdomains of dimensions $`20\times5\times5`$ with different material properties. Both sides respect a same constitutive law and a mapping is applied on the material properties.

The beam is clamped on its left face (Dirichlet boundary condition) and a vertical traction force is applied on its right face (Neumann Boundary condition).

The objective is to find the resulting displacement. Even though the proposed problem is linear, a non-linear solver is used for the example. Linear solvers have been proposed for the Stokes and Thermodynamic problems.

## Implementation

### Libraries
Computing the Finite element problem within FEniCSx in python requires to load the libraries:

```python
import dolfinx
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import ufl
import basix
import petsc4py
import mpi4py
import numpy
import pyvista
``` 
One can assess the version of FEniCSx with the following:
```python
print("Dolfinx version is:",dolfinx.__version__)
```
### Mesh generation

FEniCSx allows the creation of rectangle and boxes directly within its framework. It is however recommended to use **GMSH** to generate more complex structures since it exists a strong compatibility between GMSH and FEniCSx.

```python
domain = dolfinx.mesh.create_box(mpi4py.MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [40, 1, 1]], [20, 5, 5], dolfinx.mesh.CellType.hexahedron)
```

Once the mesh is defined, we can identify the subdomains. Locators (functions of space) and markers (tags) need to be introduced. For instance, to identify the cells from the subdomains, we define the following locators:

```python
def Omega_left(x):
    return x[0] <= 0.5*L
# 
def Omega_right(x):
    return x[0] >= 0.5*L
```
Once the locators are defined, we can identify the indices of the cells based on their position:
```python
cells_left  = dolfinx.mesh.locate_entities(domain, domain.topology.dim, Omega_left)
cells_right = dolfinx.mesh.locate_entities(domain, domain.topology.dim, Omega_right)
```

The identification of the boundaries follows exactly the same concept. Using `locate_entities_boundary` allows to create the connectivity.

```python
# Boundary locators
def left(x):
    return numpy.isclose(x[0], 0)
# 
def right(x):
    return numpy.isclose(x[0], L)
# 
def bottom(x):
    return numpy.isclose(x[2], 0)
# 
# Mark the boundaries
fdim          = domain.topology.dim - 1
left_facets   = dolfinx.mesh.locate_entities_boundary(domain, fdim, left)
right_facets  = dolfinx.mesh.locate_entities_boundary(domain, fdim, right)
bottom_facets = dolfinx.mesh.locate_entities_boundary(domain, fdim, bottom)
# 
# Concatenate and sort the arrays based on facet indices. Left facets marked with 1, right facets with two
marked_facets = numpy.hstack([left_facets, right_facets, bottom_facets])
marked_values = numpy.hstack([numpy.full_like(left_facets, 1), numpy.full_like(right_facets, 2), numpy.full_like(bottom_facets, 3)])
sorted_facets = numpy.argsort(marked_facets)
facet_tag     = dolfinx.mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])
```

To assess that the domain is well tagged, an XDMF file can be created as follows:

```python
with dolfinx.io.XDMFFile(mpi4py.MPI.COMM_WORLD, "tags.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_meshtags(facet_tag,domain.geometry)
```

### Material parameters
This example relies on a multimaterial definition based on the mapping of the material parameters. To do so, a DG0 function (defined at the Gauss points) attributes a Young modulus value to each cell based on its location:

```python3
DG0_space = dolfinx.fem.functionspace(domain, ("DG", 0))
# Map the Young's Modulus
E                      = dolfinx.fem.Function(DG0_space)
E.x.array[cells_left]  = numpy.full_like(cells_left, 1e8, dtype=dolfinx.default_scalar_type)
E.x.array[cells_right] = numpy.full_like(cells_right, 2.5e4, dtype=dolfinx.default_scalar_type)
```
The Poisson ratio has been kept constant for all subdomains:

```python3
nu = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type(0.3))
```

A mapping of the Lamé coefficients is then proposed by:
```python3
# Lamé Coefficients
lmbda_m        = E*nu.value/((1+nu.value)*(1-2*nu.value))   
mu_m           = E/(2*(1+nu.value)) 
```

### Function spaces, Functions and operators

To identify the displacement, we chose a vectorial 2nd order Lagrange representation (P2). The XDMF does not support high order functions so we also create a first order space in which we will interpolate the solution:
```python3
# Vector Element
P1_v = basix.ufl.element("P", domain.topology.cell_name(), degree=1, shape=(domain.topology.dim,))
P2_v = basix.ufl.element("P", domain.topology.cell_name(), degree=2, shape=(domain.topology.dim,))
# Function_spaces
P1v_space = dolfinx.fem.functionspace(domain, P1_v)
V         = dolfinx.fem.functionspace(domain, P2_v)
```

The mathematical spaces being defined, one can introduce the functions, expressions for interpolation, test functions and trial functions. It is recommended to place them all at a same position for debugging.

```python3
v  = ufl.TestFunction(V)
u  = dolfinx.fem.Function(V)
du = ufl.TrialFunction(V)
u_export      = dolfinx.fem.Function(P1v_space)
u_export.name = "u"
u_expr        = dolfinx.fem.Expression(u,P1v_space.element.interpolation_points())
u_export.interpolate(u_expr)
u_export.x.scatter_forward()
```
To evaluate a reaction force or a displacement over a surface, a form can be used such that:
```python3
# Evaluation of the displacement on the edge
Nz                = dolfinx.fem.Constant(domain, numpy.asarray((0.0,0.0,1.0)))
Displacement_expr = dolfinx.fem.form((ufl.dot(u,Nz))*ds(2))
```
is equivalent to:
```math
\frac{1}{S}\int u\cdot Nz \mathrm{d}S
```

For a volume, we would have had $`\frac{1}{V}\int f \mathrm{d}\Omega`$ computed with:
```python3
volume_eval = dolfinx.fem.form(f*dx)
```

The form is computed later after the solver application. 