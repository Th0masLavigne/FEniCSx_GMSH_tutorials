import dolfinx
import mpi4py
import numpy
import pyvista
import ufl
import basix

# from dolfinx import fem, mesh, plot
L = 20.0
domain = dolfinx.mesh.create_box(mpi4py.MPI.COMM_WORLD, [[0.0, 0.0, 0.0], [L, 1, 1]], [20, 5, 5], dolfinx.mesh.CellType.hexahedron)
V = dolfinx.fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim, )))

def left(x):
    return numpy.isclose(x[0], 0)

# 
def right(x):
    return numpy.isclose(x[0], L)
# 
def bottom(x):
    return numpy.isclose(x[2], 0)
# 
fdim = domain.topology.dim - 1
left_facets = dolfinx.mesh.locate_entities_boundary(domain, fdim, left)
right_facets = dolfinx.mesh.locate_entities_boundary(domain, fdim, right)
bottom_facets = dolfinx.mesh.locate_entities_boundary(domain, fdim, bottom)

# Concatenate and sort the arrays based on facet indices. Left facets marked with 1, right facets with two
marked_facets = numpy.hstack([left_facets, right_facets, bottom_facets])
marked_values = numpy.hstack([numpy.full_like(left_facets, 1), numpy.full_like(right_facets, 2), numpy.full_like(bottom_facets, 3)])
sorted_facets = numpy.argsort(marked_facets)
facet_tag = dolfinx.mesh.meshtags(domain, fdim, marked_facets[sorted_facets], marked_values[sorted_facets])
# 

with dolfinx.io.XDMFFile(mpi4py.MPI.COMM_WORLD, "verif.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_meshtags(facet_tag,domain.geometry)



# 
u_bc = numpy.array((0,) * domain.geometry.dim, dtype=dolfinx.default_scalar_type)
# 
left_dofs = dolfinx.fem.locate_dofs_topological(V, facet_tag.dim, facet_tag.find(1))
bcs = [dolfinx.fem.dirichletbc(u_bc, left_dofs, V)]

B = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
T = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
Penalty = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type(100))

v = ufl.TestFunction(V)
u = dolfinx.fem.Function(V)
du = ufl.TrialFunction(V)

# Vector Element
P1_v     = basix.ufl.element("P", domain.topology.cell_name(), degree=1, shape=(domain.topology.dim,))
P1v_space = dolfinx.fem.functionspace(domain, P1_v)
u_export      = dolfinx.fem.Function(P1v_space)
u_export.name = "u"

# Spatial dimension
d = len(u)

# Identity tensor
I = ufl.variable(ufl.Identity(d))

# Deformation gradient
F = ufl.variable(I + ufl.grad(u))

# Right Cauchy-Green tensor
C = ufl.variable(F.T * F)

# Invariants of deformation tensors
Ic = ufl.variable(ufl.tr(C))
J = ufl.variable(ufl.det(F))

# Elasticity parameters
E = dolfinx.default_scalar_type(1.0e4)
nu = dolfinx.default_scalar_type(0.3)
mu = dolfinx.fem.Constant(domain, E / (2 * (1 + nu)))
lmbda = dolfinx.fem.Constant(domain, E * nu / ((1 + nu) * (1 - 2 * nu)))
# Stored strain energy density (compressible neo-Hookean model)
psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J))**2
# Stress
# Hyper-elasticity
P = ufl.diff(psi, F)

# P = 2.0 * mu * ufl.sym(ufl.grad(u)) + lmbda * ufl.tr(ufl.sym(ufl.grad(u))) * I

metadata = {"quadrature_degree": 4}
ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tag, metadata=metadata)
dx = ufl.Measure("dx", domain=domain, metadata=metadata)
# Define form F (we want to find u such that F(u) = 0)
Form = ufl.inner(ufl.grad(v), P) * dx - ufl.inner(v, B) * dx - ufl.inner(v, T) * ds(2)
# Contact
Form += Penalty*ufl.dot(v[2],(u[2]-(-4))) * ds(3)
# F += ufl.conditional((u[2]<-4),Penalty,0)*ufl.dot(v[2],(u[2]-(-4))) * ds(3)
Jd = ufl.derivative(Form, u, du)

from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import petsc4py
# 
problem = NonlinearProblem(Form, u, bcs, J=Jd)
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
# 
ksp = solver.krylov_solver
opts = petsc4py.PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
# Plus rapide
# opts[f"{option_prefix}ksp_type"] = "preonly"
# opts[f"{option_prefix}pc_type"] = "lu"
# opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
# Plus lent
opts[f"{option_prefix}ksp_type"] = "cg"
# opts[f"{option_prefix}pc_type"] = "gamg"
opts[f"{option_prefix}pc_type"] = "lu"
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
ksp.setFromOptions()


pyvista.start_xvfb()
plotter = pyvista.Plotter()
plotter.open_gif("deformation.gif", fps=3)

topology, cells, geometry = dolfinx.plot.vtk_mesh(u.function_space)
function_grid = pyvista.UnstructuredGrid(topology, cells, geometry)

values = numpy.zeros((geometry.shape[0], 3))
values[:, :len(u)] = u.x.array.reshape(geometry.shape[0], len(u))
function_grid["u"] = values
function_grid.set_active_vectors("u")

# Warp mesh by deformation
warped = function_grid.warp_by_vector("u", factor=1)
warped.set_active_vectors("u")

# Add mesh to plotter and visualize
actor = plotter.add_mesh(warped, show_edges=True, lighting=False, clim=[0, 10])

# Compute magnitude of displacement to visualize in GIF
Vs = dolfinx.fem.functionspace(domain, ("Lagrange", 2))
magnitude = dolfinx.fem.Function(Vs)
us = dolfinx.fem.Expression(ufl.sqrt(sum([u[i]**2 for i in range(len(u))])), Vs.element.interpolation_points())
magnitude.interpolate(us)
warped["mag"] = magnitude.x.array

from dolfinx import log
log.set_log_level(log.LogLevel.INFO)


u_expr = dolfinx.fem.Expression(u,P1v_space.element.interpolation_points())
u_export.interpolate(u_expr)
u_export.x.scatter_forward()
t=0
tval0 = -0.75
# tval0 = -0.5
# for n in range(1, 10):
for n in range(1, 10):
    T.value[2] = n * tval0
    num_its, converged = solver.solve(u)
    assert (converged)
    u.x.scatter_forward()
    print(min(u.x.array[:]))
    print(f"Time step {n}, Number of iterations {num_its}, Load {T.value}")
    function_grid["u"][:, :len(u)] = u.x.array.reshape(geometry.shape[0], len(u))
    magnitude.interpolate(us)
    warped.set_active_scalars("mag")
    warped_n = function_grid.warp_by_vector(factor=1)
    warped.points[:, :] = warped_n.points
    warped.point_data["mag"][:] = magnitude.x.array
    plotter.update_scalar_bar_range([0, 10])
    plotter.write_frame()
    xdmf.write_function(u_export,t+tval0)
plotter.close()
xdmf.close()

