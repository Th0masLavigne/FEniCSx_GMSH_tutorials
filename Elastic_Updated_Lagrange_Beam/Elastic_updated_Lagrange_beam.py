# Thomas Lavigne
# 02-08-2024
# 
#----------------------------------------------------------------------
# Libraries
#----------------------------------------------------------------------
# 
import dolfinx
import mpi4py
import numpy
import pyvista
import ufl
import basix
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import petsc4py
# 
print("Dolfinx version is:",dolfinx.__version__)
# 
#----------------------------------------------------------------------
# Loading of the FE mesh
#----------------------------------------------------------------------
# 
# Create the mesh in dolfinx and identify the volumes / Boundaries
L      = 40.0
domain = dolfinx.mesh.create_box(mpi4py.MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [L, 1, 1]], [20, 5, 5], dolfinx.mesh.CellType.hexahedron)
# 
# 
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
# 
# 
with dolfinx.io.XDMFFile(mpi4py.MPI.COMM_WORLD, "tags.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_meshtags(facet_tag,domain.geometry)
# 
# 
#----------------------------------------------------------------------
# Problem and Material data
#----------------------------------------------------------------------
# Map the Young's Modulus
E = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type(1e3))
# 
# Poisson ratio
nu = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type(0.1))
# 
# Lamé Coefficients
lmbda_m        = E*nu.value/((1+nu.value)*(1-2*nu.value))   
mu_m           = E/(2*(1+nu.value)) 
# 
#----------------------------------------------------------------------
# Definition of functional spaces
#----------------------------------------------------------------------
#
# Vector Element
P1_v = basix.ufl.element("P", domain.topology.cell_name(), degree=1, shape=(domain.topology.dim,))
P2_v = basix.ufl.element("P", domain.topology.cell_name(), degree=2, shape=(domain.topology.dim,))
# Function_spaces
P1v_space = dolfinx.fem.functionspace(domain, P1_v)
V         = dolfinx.fem.functionspace(domain, P2_v)
# 
updated_mesh_space = dolfinx.fem.functionspace(domain, domain.ufl_domain().ufl_coordinate_element())
# 
#----------------------------------------------------------------------
# Functions
v  = ufl.TestFunction(V)
u  = dolfinx.fem.Function(V)
u_n  = dolfinx.fem.Function(V)
du = ufl.TrialFunction(V)
# 
u_export      = dolfinx.fem.Function(P1v_space)
du_update     = dolfinx.fem.Function(updated_mesh_space)
u_export.name = "u"
#----------------------------------------------------------------------
# Operators
metadata = {"quadrature_degree": 4}
ds       = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
dx       = ufl.Measure("dx", domain=domain, metadata=metadata)
#----------------------------------------------------------------------
# Expressions
# XDMF needs linear displacement
u_expr   = dolfinx.fem.Expression(u_n,P1v_space.element.interpolation_points())
u_export.interpolate(u_expr)
u_export.x.scatter_forward()
# 
# Evaluation of the displacement on the edge
Nx                = dolfinx.fem.Constant(domain, numpy.asarray((1.0,0.0,0.0)))
Displacement_expr = dolfinx.fem.form((ufl.dot(u,Nx))*ds(3))
# 
#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------
# 
u_bc = numpy.array((0,) * domain.geometry.dim, dtype=dolfinx.default_scalar_type)
# 
left_dofs = dolfinx.fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))
bcs       = [dolfinx.fem.dirichletbc(u_bc, left_dofs, V)]
# 
#----------------------------------------------------------------------
# Definition of the Variationnal Form
#----------------------------------------------------------------------
# 
# Body forces vector
B = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
# Traction force vector
T = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
# 
# Constitutive Law
def Hookean(mu,lmbda,u):
    return 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * ufl.variable(ufl.Identity(len(u)))
# 
#----------------------------------------------------------------------
# Variational form
# Define form F (we want to find u such that F(u) = 0)
F   = ufl.inner(ufl.grad(v), Hookean(mu_m,lmbda_m,u_n+u)) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)
# 
# Jacobian of the problem
J__ = ufl.derivative(F, u, du)
# 
#----------------------------------------------------------------------
# Problem and Solver settings
problem = NonlinearProblem(F, u, bcs, J=J__)
solver = NewtonSolver(domain.comm, problem)
# 
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
# 
# 
#----------------------------------------------------------------------
# Solving and post-processing
#----------------------------------------------------------------------
# 
#----------------------------------------------------------------------
# Post-processing
pyvista.start_xvfb()
plotter = pyvista.Plotter()
plotter.open_gif("elastic_UL_beam.gif", fps=3)
# 
topology, cells, geometry = dolfinx.plot.vtk_mesh(u.function_space)
function_grid             = pyvista.UnstructuredGrid(topology, cells, geometry)
# 
values             = numpy.zeros((geometry.shape[0], 3))
values[:, :len(u)] = u_n.x.array.reshape(geometry.shape[0], len(u_n))
function_grid["u"] = values
function_grid.set_active_vectors("u")
# 
# Warp mesh by deformation
warped        = function_grid.warp_by_vector("u", factor=1)
warped.set_active_vectors("u")
# 
# Add mesh to plotter and visualize
actor         = plotter.add_mesh(warped, show_edges=True, lighting=False, clim=[0, 10])
# 
# Compute magnitude of displacement to visualize in GIF
Vs            = dolfinx.fem.functionspace(domain, ("Lagrange", 2))
magnitude     = dolfinx.fem.Function(Vs)
us            = dolfinx.fem.Expression(ufl.sqrt(sum([u_n[i]**2 for i in range(len(u_n))])), Vs.element.interpolation_points())
magnitude.interpolate(us)
warped["mag"] = magnitude.x.array
# 
#----------------------------------------------------------------------
# Debug instance
log_solve=True
if log_solve:
    from dolfinx import log
    log.set_log_level(log.LogLevel.INFO)
#----------------------------------------------------------------------
# Computation (an increasing load allows to update the initial condition)
# Load increment
tval0 = 1
# Loop to get to the total load
for n in range(1, 10):
    T.value[0] = n*tval0
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
    u_n.x.array[:]+=u.x.array[:]
    u_n.x.scatter_forward()
    du_update.interpolate(u)
    du_update.x.scatter_forward()
    # Evaluate the displacement
    displacement_      = dolfinx.fem.assemble_scalar(Displacement_expr)
    Surface            = 1*1
    displacement_right = 1/Surface*domain.comm.allreduce(displacement_, op=mpi4py.MPI.SUM)
    print("Edge displacement increment:", displacement_right)
    # 
    print(f"Time step {n}, Number of iterations {num_its}, Load {T.value}")
    # Post-processing
    function_grid["u"][:, :len(u_n)] = u_n.x.array.reshape(geometry.shape[0], len(u_n))
    magnitude.interpolate(us)
    warped.set_active_scalars("mag")
    warped_n                    = function_grid.warp_by_vector(factor=10)
    warped.points[:, :]         = warped_n.points
    warped.point_data["mag"][:] = magnitude.x.array
    plotter.update_scalar_bar_range([0, 0.05])
    plotter.write_frame()
    xdmf.write_function(u_export,n*tval0)
    domain.geometry.x[:, :domain.geometry.dim] += du_update.x.array.reshape((-1, domain.geometry.dim))
plotter.close()
xdmf.close()
# EoF