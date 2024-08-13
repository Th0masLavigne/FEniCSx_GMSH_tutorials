from fenics import *
from mshr import *
from dolfin import *


#----------------------------------------------------------------------
# Problem and Material data
#----------------------------------------------------------------------
T_i   = 100.0
T_e   = 20.0
kc    = 1.0
hc    = 10.0
f     = 0.
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Geometrical parameters
#----------------------------------------------------------------------
r_c = 0.02
l_d = 0.24
h_d = 0.1
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Construction of the FE mesh
#----------------------------------------------------------------------
C1     = Circle(Point(0.,0.), r_c)
R1     = Rectangle(Point(-0.5*l_d, 0.), Point(0.5*l_d, h_d))
domain = R1 - C1
mesh   = generate_mesh(domain, 50)

# Mesh writting in file
file_mesh = File('mesh.pvd')
file_mesh << mesh
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of relefant mesh boundaries
#----------------------------------------------------------------------
bound_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
bound_markers.set_all(0)
tol = 1.0e-3

# Upper bound (NB: quando definisco una classe lasciare una linea)
class Boundary_up(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], h_d, tol)

bound_up = Boundary_up()
bound_up.mark(bound_markers, 1)

# Circle bound
class Boundary_c(SubDomain):
	def inside(self, x, on_boundary):
		return abs(sqrt(x[1]*x[1]+x[0]*x[0])-r_c) < tol and on_boundary

bound_c = Boundary_c()
bound_c.mark(bound_markers, 2)

# Redefiniction dx and ds
dx = Measure("dx", domain = mesh)
ds = Measure("ds", domain = mesh, subdomain_data = bound_markers)

# Boundaries writting in file
file_boundary = File('boundary.pvd')
file_boundary << bound_markers
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of functional spaces
#----------------------------------------------------------------------
V = FunctionSpace(mesh, 'P', 1)
u = TrialFunction(V)
v = TestFunction(V)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of Dirichlet boundary conditions
# dbc_1 = DirichletBC(V, valeur, ligne, marker)
# dbc   = [dbc_1, dbc_2, ..., dbc_n]
#----------------------------------------------------------------------
dbc_1 = DirichletBC(V, T_i, bound_markers, 2)
dbc   = [dbc_1]
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of the Variationnal Form
#----------------------------------------------------------------------
#A = kc*dot(grad(u), grad(v))*dx 
#L = f*v*dx 
#----------------------------------------------------------------------
A = kc*dot(grad(u), grad(v))*dx + hc*u*v*ds(1)
L = f*v*dx + hc*T_e*v*ds(1)

# Alternative method
#F= kc*dot(grad(u), grad(v))*dx - f*v*dx +hc*(u-T_e)*v*ds(1)
#A,L=lhs(F),rhs(F)
#----------------------------------------------------------------------


# Solution
T = Function(V)
solve(A == L, T, dbc)

file_T = File('Temperature.pvd')
file_T << T




