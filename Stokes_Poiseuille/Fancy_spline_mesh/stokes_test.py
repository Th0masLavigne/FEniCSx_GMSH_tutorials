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
from dolfinx.fem.petsc import LinearProblem
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
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh('./2D_Stokes_fancy_mesh.msh', mpi4py.MPI.COMM_WORLD, 0, gdim=2)
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
sol    = ufl.TrialFunction(CHS)
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
ds    = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
# 
#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------
# 
bcs = []
fdim = mesh.topology.dim - 1
# 
def add_dirichlet_BC(functionspace,dimension,facet,value):
    dofs   = dolfinx.fem.locate_dofs_topological(functionspace, dimension, facet)
    bcs.append(dolfinx.fem.dirichletbc(value, dofs, functionspace))
# 
# No-slip boundary condition for velocity
# bottom
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(4), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(4), petsc4py.PETSc.ScalarType(0.))
# top
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(2), petsc4py.PETSc.ScalarType(0.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(2), petsc4py.PETSc.ScalarType(0.))
# 
# Inflow boundary condition for velocity
# left
add_dirichlet_BC(CHS.sub(1).sub(0),fdim,facet_tag.find(1), petsc4py.PETSc.ScalarType(1.))
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(1), petsc4py.PETSc.ScalarType(0.))
# 
# Outflow vy = 0
add_dirichlet_BC(CHS.sub(1).sub(1),fdim,facet_tag.find(3), petsc4py.PETSc.ScalarType(0.))
# 
#----------------------------------------------------------------------
# Weak form
#----------------------------------------------------------------------
# Define variational problem
# 
# 
# 
# Alternative form (from the tutorial to comprae with the iterative solution)
A1 = ufl.inner(ufl.grad(u), ufl.grad(w))*dx + ufl.div(w)*p*dx 
A2 = q*ufl.div(u)*dx
# Assembling of the system of eqs 
A = A1 + A2
# 
f = dolfinx.fem.Constant(mesh,(0.0, 0.0))
L = ufl.inner(f, w)*dx
# 
#----------------------------------------------------------------------
# Debug instance
log_solve=True
if log_solve:
	from dolfinx import log
	log.set_log_level(log.LogLevel.INFO)
# 
#----------------------------------------------------------------------
# Solve
#----------------------------------------------------------------------
#
# 
problem = LinearProblem(A, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()
# 
#----------------------------------------------------------------------
# Post-process
#----------------------------------------------------------------------
#
# Get sub-functions
p_, u_ = uh.split()
p_.name = "p"
# 
u_expr = dolfinx.fem.Expression(uh.sub(1),P1v_space.element.interpolation_points())
u_export.interpolate(u_expr)
u_export.x.scatter_forward()
# 

# 
xdmf = dolfinx.io.XDMFFile(mesh.comm, "2D_Stokes.xdmf", "w")
xdmf.write_mesh(mesh)
t=0
xdmf.write_function(u_export,t)
xdmf.write_function(p_,t)
xdmf.close()
# 
# EoF