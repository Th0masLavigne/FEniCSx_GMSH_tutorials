#----------------------------------------------------------------------
#---------------------------------------------------------------------- 
# Stokes equations
#----------------------------------------------------------------------
#----------------------------------------------------------------------
from fenics import *
from mshr   import *
from dolfin import *


#----------------------------------------------------------------------
# Material parameters 
#----------------------------------------------------------------------
mu  = 0.001
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Geometrical parameters
#----------------------------------------------------------------------
l_d = 4.0
h_d = 1.0
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Construction of the FE mesh
#----------------------------------------------------------------------
R1     = Rectangle(Point(0., 0.), Point(l_d, h_d))
domain = R1
mesh   = generate_mesh(domain, 80)

# Mesh writting in file
file_mesh = File('mesh.pvd')
file_mesh << mesh
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Build function space (vectorial and scalar)
#----------------------------------------------------------------------
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of relevant bounds
# This time with boundary markers to extract subdomains
#----------------------------------------------------------------------
boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
boundary_markers.set_all(0)
tol = 1.0e-3


# Right bound
class Boundary_right(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], l_d, tol)

bound_right = Boundary_right()
bound_right.mark(boundary_markers, 2)

# Left bound
class Boundary_left(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], 0, tol)

bound_left = Boundary_left()
bound_left.mark(boundary_markers, 4)

# Botton bound
class Boundary_right(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0, tol)

bound_right = Boundary_right()
bound_right.mark(boundary_markers, 1)

# Top bound
class Boundary_left(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], h_d, tol)

bound_left = Boundary_left()
bound_left.mark(boundary_markers, 3)


# Boundaries writting in file
file_boundary = File('boundary.pvd')
file_boundary << boundary_markers

# Redefiniction dx and ds
dx = Measure("dx", domain = mesh)
ds = Measure("ds", domain = mesh, subdomain_data = boundary_markers)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------

v_null = Constant((0.0, 0.0))
zero   = Constant(0)
one    = Constant(1)

# No-slip boundary condition for velocity
bc0    = DirichletBC(W.sub(0), v_null, boundary_markers, 1)
bc1    = DirichletBC(W.sub(0), v_null, boundary_markers, 3)

# Inflow boundary condition for velocity
v_prof = Expression(("sin(x[1]*pi)"), degree=2)
#bc2 = DirichletBC(W.sub(0).sub(0), v_prof, boundary_markers, 4)
bc2 = DirichletBC(W.sub(0).sub(0), one, boundary_markers, 4)
bc3 = DirichletBC(W.sub(0).sub(1), zero,   boundary_markers, 4)

# Outflow vy = 0
bc4 = DirichletBC(W.sub(0).sub(1), zero, boundary_markers, 2)

# Collect boundary conditions
bcs = [bc0, bc1, bc2, bc3, bc4]
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Definition of the Variationnal Form
#----------------------------------------------------------------------
# The bilinear and linear forms corresponding to the weak mixed
# formulation of the Stokes equations are defined as follows:
#----------------------------------------------------------------------
# Define variational problem
(u, p) = TrialFunctions(W)
(w, q) = TestFunctions(W)
b      = Constant((0.0, 0.0))
n      = FacetNormal(mesh)
p0     = 10

import dolfin as df
I = df.Identity(2)

# Divergence form
a1   = inner((- p*I + mu*grad(u) + mu*grad(u).T), grad(w))*dx
a2   = div(u)*q*dx
a    = a1 + a2

L  = dot(b, w)*dx - p0*dot(n, w)*ds(2)

# Assemble system
# A, bb = assemble_system(a, L, bcs)
U = Function(W)
solve(a == L, U, bcs, solver_parameters={'linear_solver':'lu'})


# Get sub-functions
u, p = U.split()

# Save solution in VTK format
ufile_pvd = File("velocity.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << p




