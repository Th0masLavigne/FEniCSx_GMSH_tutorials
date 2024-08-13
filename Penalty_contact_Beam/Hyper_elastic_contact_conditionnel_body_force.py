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
L      = 20.0
domain = dolfinx.mesh.create_box(mpi4py.MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [L, 1, 1]], [20, 5, 5], dolfinx.mesh.CellType.hexahedron)
# 
# Locators
# 
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
#  
# Elasticity parameters
E     = dolfinx.default_scalar_type(1.0e4)
nu    = dolfinx.default_scalar_type(0.3)
mu    = dolfinx.fem.Constant(domain, E / (2 * (1 + nu)))
lmbda = dolfinx.fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
Penalty = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type(1000))
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
#----------------------------------------------------------------------
# Functions
v  = ufl.TestFunction(V)
u  = dolfinx.fem.Function(V)
du = ufl.TrialFunction(V)
#----------------------------------------------------------------------
# Operators
metadata = {"quadrature_degree": 4}
ds       = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
dx       = ufl.Measure("dx", domain=domain, metadata=metadata)
#----------------------------------------------------------------------
# Expressions
# XDMF needs linear displacement
u_export      = dolfinx.fem.Function(P1v_space)
u_export.name = "u"
u_expr        = dolfinx.fem.Expression(u,P1v_space.element.interpolation_points())
u_export.interpolate(u_expr)
u_export.x.scatter_forward()
# Evaluation of the displacement on the edge
Nz                = dolfinx.fem.Constant(domain, numpy.asarray((0.0,0.0,1.0)))
Displacement_expr = dolfinx.fem.form((ufl.dot(u,Nz))*ds(2))
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
# Constitutive Laws
def Neo_Hoolean(mu,lmbda):
    # Spatial dimension
    d   = len(u)
    # Identity tensor
    I   = ufl.variable(ufl.Identity(d))
    # Deformation gradient
    F   = ufl.variable(I + ufl.grad(u))
    # Right Cauchy-Green tensor
    C   = ufl.variable(F.T * F)
    # Invariants of deformation tensors
    Ic  = ufl.variable(ufl.tr(C))
    J   = ufl.variable(ufl.det(F))
    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
    return ufl.diff(psi, F)
# 
#----------------------------------------------------------------------
# Variational form
# Define form F (we want to find u such that F(u) = 0)
Form = ufl.inner(ufl.grad(v), Neo_Hoolean(mu,lmbda)) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)
# Contact based on the gap if the gap is negative, 0 otherwise
Form += Penalty*ufl.dot(v[2],ufl.conditional((u[2]<-4),(u[2]-(-4)),0)) * ds(3)
# 
# Jacobian of the problem
Jd = ufl.derivative(Form, u, du)
# 
#----------------------------------------------------------------------
# Problem and Solver settings
problem = NonlinearProblem(Form, u, bcs, J=Jd)
solver  = NewtonSolver(domain.comm, problem)
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
plotter.open_gif("hyper_elastic_penalty_beam.gif", fps=10)
# 
topology, cells, geometry = dolfinx.plot.vtk_mesh(u.function_space)
function_grid             = pyvista.UnstructuredGrid(topology, cells, geometry)
# 
values             = numpy.zeros((geometry.shape[0], 3))
values[:, :len(u)] = u.x.array.reshape(geometry.shape[0], len(u))
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
us            = dolfinx.fem.Expression(ufl.sqrt(sum([u[i]**2 for i in range(len(u))])), Vs.element.interpolation_points())
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
bval0 = -0.15
# Loop to get to the total load
for n in range(1, 80):
    B.value[2]         = n * bval0
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
    print(f"Time step {n}, Number of iterations {num_its}, Load {B.value}")
    # Post-processing
    function_grid["u"][:, :len(u)] = u.x.array.reshape(geometry.shape[0], len(u))
    magnitude.interpolate(us)
    warped.set_active_scalars("mag")
    warped_n                    = function_grid.warp_by_vector(factor=1)
    warped.points[:, :]         = warped_n.points
    warped.point_data["mag"][:] = magnitude.x.array
    plotter.update_scalar_bar_range([0, numpy.max(magnitude.x.array[:])])
    plotter.write_frame()
    xdmf.write_function(u_export,n*bval0)
plotter.close()
xdmf.close()
# EoF