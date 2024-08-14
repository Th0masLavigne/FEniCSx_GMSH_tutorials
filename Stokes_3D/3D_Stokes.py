# Thomas Lavigne
# 02-08-2024
# 
#----------------------------------------------------------------------
# Libraries
#----------------------------------------------------------------------
# 
import dolfinx
import ufl
import basix
import mpi4py
import petsc4py
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
# 
print("Dolfinx version is:",dolfinx.__version__)
# 
#----------------------------------------------------------------------
# Loading of the FE mesh
#----------------------------------------------------------------------
# 
# Alternative for tags identification : see https://doi.org/10.1016/j.jmbbm.2023.105902 (v 0.5.2) and https://jsdokken.com/dolfinx-tutorial/
# The alternative method is useful in case of a mesh which is not already tagged.
# 
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh('./3D_Stokes_mesh.msh', mpi4py.MPI.COMM_WORLD, 0, gdim=3)
# 
#----------------------------------------------------------------------
# Definition of functional space, trial & test functions
#----------------------------------------------------------------------
#___________________________________________________________________________
# Define function space
# Finite Element 
P1       = basix.ufl.element("P", mesh.topology.cell_name(), degree=1)
# Vector Element
P1_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=1, shape=(mesh.topology.dim,))
P2_v     = basix.ufl.element("P", mesh.topology.cell_name(), degree=2, shape=(mesh.topology.dim,))
# Mixed element
MxE      = basix.ufl.mixed_element([P1,P2_v])
# 
P1v_space = dolfinx.fem.functionspace(mesh, P1_v)
CHS       = dolfinx.fem.functionspace(mesh, MxE)
# 
tensor_elem  = basix.ufl.element("P", mesh.topology.cell_name(), degree=1, shape=(3,3))
tensor_space = dolfinx.fem.functionspace(mesh, tensor_elem)
#___________________________________________________________________________
# Define function & parameters
# 
u_export      = dolfinx.fem.Function(P1v_space)
u_export.name = "u"
# 
sol    = dolfinx.fem.Function(CHS)
q, w   = ufl.TestFunctions(CHS)
# Solution vector
p, u   = ufl.split(sol)
# 
# Definition of the normal vector
n      = ufl.FacetNormal(mesh)
# 
# Specify the desired quadrature degree
q_deg = 4
# Redefinition dx and ds
dx    = ufl.Measure('dx', metadata={"quadrature_degree":q_deg}, subdomain_data=cell_tag, domain=mesh)
ds    = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag, metadata={"quadrature_degree":q_deg})
# 
#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------
# 
# tag_left, tag_top_left, tag_middle_up, tag_top_right, tag_right, tag_bottom, tag front = 1, 2, 3, 4, 5, 6, 7
# 
bcs = []
fdim = mesh.topology.dim - 1
# 
def add_dirichlet_BC(functionspace,dimension,facet,value):
	dofs   = dolfinx.fem.locate_dofs_topological(functionspace, dimension, facet)
	bcs.append(dolfinx.fem.dirichletbc(value, dofs, functionspace))
# 
# No-slip boundary condition for velocity
# top left
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(2), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(2), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(2),fdim,facet_tag.find(2), petsc4py.PETSc.ScalarType(0.))
# top right
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(4), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(4), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(2),fdim,facet_tag.find(4), petsc4py.PETSc.ScalarType(0.))
# middle up
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(3), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(3), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(2),fdim,facet_tag.find(3), petsc4py.PETSc.ScalarType(0.))
# 
# inflow vx = 1 left side
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(1), petsc4py.PETSc.ScalarType(1.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(1), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(2),fdim,facet_tag.find(1), petsc4py.PETSc.ScalarType(0.))
# 
# Symmetry plan vy = 0 bottom
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(6), petsc4py.PETSc.ScalarType(0.))
# 
# Symmetry plan vz = 0 front
add_dirichlet_BC(CHS.sub(1).sub(2),fdim,facet_tag.find(7), petsc4py.PETSc.ScalarType(0.))
# 
#----------------------------------------------------------------------
# Weak form
#----------------------------------------------------------------------
# Define variational problem
# 
b = dolfinx.fem.Constant(mesh,(0.0, 0.0, 0.0))
F = ufl.inner(ufl.grad(u), ufl.grad(w))*dx + ufl.div(w)*p*dx + q*ufl.div(u)*dx - ufl.inner(b, w)*dx
# 
#----------------------------------------------------------------------
# Debug instance
log_solve=True
if log_solve:
	from dolfinx import log
	log.set_log_level(log.LogLevel.INFO)
# 
#----------------------------------------------------------------------
problem = NonlinearProblem(F, sol, bcs)
solver = NewtonSolver(mesh.comm, problem)
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
num_its, converged = solver.solve(sol)
sol.x.scatter_forward()
# Get sub-functions
p_, u_ = sol.split()
p_.name = "p"
# 
u_expr = dolfinx.fem.Expression(sol.sub(1),P1v_space.element.interpolation_points())
u_export.interpolate(u_expr)
u_export.x.scatter_forward()
# 
strainrate=dolfinx.fem.Function(tensor_space)
strainrate.name = "strainrate"
# 0.5*(grad_cyl(u_export) + grad_cyl(u_export).T)
strainrate_expr = dolfinx.fem.Expression(ufl.sym(ufl.grad(u_)),tensor_space.element.interpolation_points())
strainrate.interpolate(strainrate_expr)
strainrate.x.scatter_forward()
# # 
Id = ufl.Identity(3)
eta=1.
stress=dolfinx.fem.Function(tensor_space)
stress.name = "stress"
stress_expr = dolfinx.fem.Expression(-1.*p_*Id + eta*2*ufl.sym(ufl.grad(u_)),tensor_space.element.interpolation_points())
stress.interpolate(stress_expr)
stress.x.scatter_forward()
# 
xdmf = dolfinx.io.XDMFFile(mesh.comm, "3D_Stokes.xdmf", "w")
xdmf.write_mesh(mesh)
t=0
xdmf.write_function(u_export,t)
xdmf.write_function(p_,t)
xdmf.write_function(strainrate,t)
xdmf.write_function(stress,t)
# 
import pyvista
import numpy
pyvista.start_xvfb()
topology, cell_types, geometry = dolfinx.plot.vtk_mesh(P1v_space)
values = numpy.zeros((geometry.shape[0], 3), dtype=numpy.float64)
values[:, :len(u_export)] = u_export.x.array.real.reshape((geometry.shape[0], len(u_export)))

# Create a point cloud of glyphs
function_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
function_grid["u"] = values
glyphs = function_grid.glyph(orient="u", factor=0.1)

# Create a pyvista-grid for the mesh
grid = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(mesh, mesh.topology.dim))

# Create plotter
plotter = pyvista.Plotter()
# plotter.add_mesh(grid, style="wireframe", color="k")
plotter.add_mesh(grid, color="grey", ambient=0.02, opacity=0.05 , show_edges=False)
plotter.add_mesh(glyphs, cmap='coolwarm')
plotter.update_scalar_bar_range([1, 8])
plotter.view_vector((0, 3, -20))
plotter.show_axes()
# plotter.view_xy()
plotter.save_graphic('Result.pdf')
plotter.close()
# EoF