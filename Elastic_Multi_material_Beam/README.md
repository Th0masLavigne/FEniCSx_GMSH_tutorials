# Elastic Multimaterial Beam

## Description of the problem

This example aims in providing an example of a 3D beam of dimensions $`40\times1\times1`$ composed of two materials.

The beam is subdivided into two subdomains of dimensions $`20\times1\times1`$ with different material properties. Both sides respect a same constitutive law and a mapping is applied on the material properties.

![Geometry](./Elastic_multimaterial_beam.jpg)

The beam is clamped on its left face (Dirichlet boundary condition) and a vertical traction force is applied on its right face (Neumann Boundary condition).

The objective is to find the resulting displacement. 

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
L = 40
domain = dolfinx.mesh.create_box(mpi4py.MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [L, 1, 1]], [20, 5, 5], dolfinx.mesh.CellType.hexahedron)
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

The identification and marking of the boundaries follows exactly the same concept. Using `locate_entities_boundary` allows to create the connectivity.

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

To verify if the domain is well tagged, an XDMF file can be created as follows:

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

The solid is assumed to follow the Hookean constitutive law such that $` \mathbf{\sigma}(\mathbf{u}) = 2 \mu \mathbf{\varepsilon}(\mathbf{u}) + \lambda \mathrm{tr}(\mathbf{\varepsilon}(\mathbf{u}))\mathbf{I_d}:`$
```python
# Constitutive Law
def Hookean(mu,lmbda):
    return 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * ufl.variable(ufl.Identity(len(u)))
```
**Remark:** Note that the hereabove function must be introduced after the definition of u.

The body forces and traction forces are defined using:
```python
# Body forces vector
B = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
# Traction force vector
T = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
```
The choice of a constant allows to dynamically update the value with time. It is of interest for boundary conditions and loading.

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

The following operators are also defined:
```python
metadata = {"quadrature_degree": 4}
ds       = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
dx       = ufl.Measure("dx", domain=domain, metadata=metadata)
```

To evaluate a reaction force or a displacement over a surface, a form can be used such that:
```python3
# Evaluation of the displacement on the edge
Nz                = dolfinx.fem.Constant(domain, numpy.asarray((0.0,0.0,1.0)))
Displacement_expr = dolfinx.fem.form((ufl.dot(u,Nz))*ds(2))
```
is equivalent to:
```math
\frac{1}{{\partial\Omega_t}}\int_{\partial\Omega_t} u\cdot N_t \mathrm{d}S,~\text{with } \overrightarrow{T}=T~N_t
```

For a volume, we would have had $`\frac{1}{V}\int f \mathrm{d}\Omega`$ computed with:
```python3
volume_eval = dolfinx.fem.form(f*dx)
```

The form is computed later after the solver application. 

### Dirichlet boundary conditions
The boundary condition being fixed (no dynamically imposed displacement), the clamp is defined as follows:

```python
u_bc = numpy.array((0,) * domain.geometry.dim, dtype=dolfinx.default_scalar_type)
# 
left_dofs = dolfinx.fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))
bcs       = [dolfinx.fem.dirichletbc(u_bc, left_dofs, V)]
```

### Variationnal form

For an elastic problem, the variationnal form to be solved is:

```math
\int_\Omega\sigma(u):\varepsilon(v)\mathrm{d}\Omega - \int_\Omega B\cdot v \mathrm{d}\Omega -  \int_{\partial\Omega_t}T\cdot v \mathrm{d}S = 0
```
where B stands for the body forces, T the traction forces, u is the unknown and v the test function. 

This is traduced in FEniCSx with:

```python
F   = ufl.inner(ufl.grad(v), Hookean(mu_m,lmbda_m)) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)
```
The Jacobian of the problem can further be defined with:

```python
J__ = ufl.derivative(F, u, du)
```

Finally the problem is introduced as:
```python
problem = NonlinearProblem(F, u, bcs, J=J__)
solver = NewtonSolver(domain.comm, problem)
```

### Solver settings:
The solver settings are defined as follows:
```python
# Absolute tolerance
solver.atol                  = 1e-8
# relative tolerance
solver.rtol                  = 1e-8
# Convergence criterion
solver.convergence_criterion = "incremental"
# Maximum iterations
solver.max_it                = 15
# Solver Pre-requisites
ksp                                               = solver.krylov_solver
opts                                              = petsc4py.PETSc.Options()
option_prefix                                     = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"]                  = "preonly"
opts[f"{option_prefix}pc_type"]                   = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()
```

### Solving and post-processing

This example provide a gif with the displacement magnitude as well as a xdmf file. Please refer to the python file. A minimal resolution code is hereafter presented.

To have the full computation log, the following is required. These information are crucial when debugging.

```python
#----------------------------------------------------------------------
# Debug instance
log_solve=True
if log_solve:
    from dolfinx import log
    log.set_log_level(log.LogLevel.INFO)
#----------------------------------------------------------------------
```

For stability concerns, we increment the load to reach the solution:

```python
# Load increment
tval0 = -0.75
# Loop to get to the total load
for n in range(1, 10):
    T.value[2]         = n * tval0
    num_its, converged = solver.solve(u)
    u.x.scatter_forward()
    try:
        assert (converged)
    except:
        if MPI.COMM_WORLD.rank == 0:
            print("*************") 
            print("Solver failed")
            print("*************") 
        break
    # 
    u_export.interpolate(u_expr)
    u_export.x.scatter_forward()
    # Evaluate the displacement
    displacement_      = dolfinx.fem.assemble_scalar(Displacement_expr)
    Surface            = 1*1
    displacement_right = 1/Surface*domain.comm.allreduce(displacement_, op=mpi4py.MPI.SUM)
    print("Edge displacement:", displacement_right)
    # 
    print(f"Time step {n}, Number of iterations {num_its}, Load {T.value}")
    xdmf.write_function(u_export,n*tval0)
xdmf.close()
```