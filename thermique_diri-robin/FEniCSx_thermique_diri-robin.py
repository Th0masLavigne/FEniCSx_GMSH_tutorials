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
mesh, cell_tag, facet_tag = dolfinx.io.gmshio.read_from_msh('./2D_thermique.msh', mpi4py.MPI.COMM_WORLD, 0, gdim=2)
# 
#----------------------------------------------------------------------
# Problem and Material data
#----------------------------------------------------------------------
T_i   = 100.0
T_e   = 20.0
kc    = 1.0
hc    = 10.0
f     = 0.
#----------------------------------------------------------------------
# 
#----------------------------------------------------------------------
# Definition of functional spaces
#----------------------------------------------------------------------
# Finite Element 
P1 = basix.ufl.element("P", mesh.topology.cell_name(), degree=1)
V  = dolfinx.fem.functionspace(mesh, P1)
u  = ufl.TrialFunction(V)
v  = ufl.TestFunction(V)
#----------------------------------------------------------------------
# 
# Specify the desired quadrature degree
q_deg = 4
# Redefinition dx and ds
dx    = ufl.Measure('dx', metadata={"quadrature_degree":q_deg}, subdomain_data=cell_tag, domain=mesh)
ds    = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tag)
# 
# 
#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------
# 
# tag_left, tag_top, tag_right, tag_bottom, tag_circle = 1, 2, 3, 4, 5
# 
bcs = []
fdim = mesh.topology.dim - 1
# 
def add_dirichlet_BC(functionspace,dimension,facet,value):
	dofs   = dolfinx.fem.locate_dofs_topological(functionspace, dimension, facet)
	bcs.append(dolfinx.fem.dirichletbc(value, dofs, functionspace))
# 
add_dirichlet_BC(V,fdim,facet_tag.find(5), petsc4py.PETSc.ScalarType(T_i))
# 
#----------------------------------------------------------------------
# Definition of the Variationnal Form
#----------------------------------------------------------------------
# 
A = kc*ufl.dot(ufl.grad(u), ufl.grad(v))*dx + hc*u*v*ds(2)
L = f*v*dx + hc*T_e*v*ds(2)
#----------------------------------------------------------------------
# Debug instance
log_solve=True
if log_solve:
	from dolfinx import log
	log.set_log_level(log.LogLevel.INFO)
# 
problem = LinearProblem(A, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()
# 
uh.name = "solution"
xdmf = dolfinx.io.XDMFFile(mesh.comm, "2D_thermique.xdmf", "w")
xdmf.write_mesh(mesh)
t=0
xdmf.write_function(uh,t)