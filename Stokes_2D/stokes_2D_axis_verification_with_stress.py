#----------------------------------------------------------------------
#---------------------------------------------------------------------- 
# Verification Stokes model axial symetry
#----------------------------------------------------------------------
#----------------------------------------------------------------------
from fenics import *
from mshr   import *
from ufl import *
from dolfin import *
from fenics import *
import random as random
import numpy as np
import time
import csv
#import matplotlib.pyplot 
from mshr import *

# Use compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True
parameters["allow_extrapolation"]           = True


#----------------------------------------------------------------------
# Model parameters
#----------------------------------------------------------------------
eta  = 1.


#----------------------------------------------------------------------
# Geometrical parameters
#----------------------------------------------------------------------
L1 = 2.
H1 = 1.
L2 = 2.
H2 = 0.5
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Construction of the FE mesh
#----------------------------------------------------------------------
R1 = Rectangle(Point(0., 0.), Point(L1, H1))
R2 = Rectangle(Point(L1, 0.), Point((L1 + L2), H2))
domain = R1 + R2;

mesh   = generate_mesh(domain, 60)

# Mesh writting in file
file_mesh = File('mesh.pvd')
file_mesh << mesh
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of functional space, trial & test functions
#----------------------------------------------------------------------
#___________________________________________________________________________
# Define function space
CP1  = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
CV2  = VectorElement("Lagrange", mesh.ufl_cell(), 2)
CHS  = FunctionSpace(mesh, MixedElement(CP1, CV2))


#___________________________________________________________________________
# Define function & parameters
sol    = TrialFunction(CHS)
q, w   = TestFunctions(CHS)
# Solution vector
p, u   = split(sol)

# Definition of the normal vector
n      = FacetNormal(mesh)
x      = SpatialCoordinate(mesh)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of relevant bounds
# This time with boundary markers to extract subdomains
#----------------------------------------------------------------------
boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
boundary_markers.set_all(0)
tol  = 1.e-2
tolc = 1.e-2

# Bounds needed for the fluid flow problem

# Boundary_1(Left bound)
class Boundary_1(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], 0, tol)

bound_1 = Boundary_1()
bound_1.mark(boundary_markers, 1)


# Boundary_2(right bound)
class Boundary_2(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], (L1+L2) , tol)

bound_2 = Boundary_2()
bound_2.mark(boundary_markers, 2)


# Boundary_3(upper bound big channel)
class Boundary_3(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], H1 , tol)

bound_3 = Boundary_3()
bound_3.mark(boundary_markers, 3)


# Boundary_4(upper bound small channel)
class Boundary_4(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], H2, tol) and x[0] >= L1

bound_4 = Boundary_4()
bound_4.mark(boundary_markers, 4)

# Boundary_5(intermediate bound)
class Boundary_5(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], L1, tol) and x[1] >= H2

bound_5 = Boundary_5()
bound_5.mark(boundary_markers, 5)


# Boundary_6 (symmetry plan)
class Boundary_6(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0., tol)

bound_6 = Boundary_6()
bound_6.mark(boundary_markers, 6)


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

v_null   = Constant((0.0, 0.0))
v_in     = Constant((1.0, 0.0))
v_prof   = Expression(("cos(x[1]*0.5*pi)"), degree=2)

zero     = Constant(0.)
p_right  = Constant(0.)

# No-slip boundary condition for velocity
bcNS3    = DirichletBC(CHS.sub(1), v_null, boundary_markers, 3)
bcNS4    = DirichletBC(CHS.sub(1), v_null, boundary_markers, 4)
bcNS5    = DirichletBC(CHS.sub(1), v_null, boundary_markers, 5)

# iuflow vx = 1
bcNS1 = DirichletBC(CHS.sub(1), v_in, boundary_markers, 1)

# Symmetry plan vy = 0
bcNS6 = DirichletBC(CHS.sub(1).sub(1), zero, boundary_markers, 6)


# Collect boundary conditions
bc_tot = [bcNS1, bcNS3, bcNS4, bcNS5, bcNS6]
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Weak form
#----------------------------------------------------------------------
# Define variational problem
import dolfin as df
I = df.Identity(3)

# Definition of divergence and gradient in cylindrycal coordinates
r = abs(x[1])

def grad_cyl(u):
	return as_tensor([[u[0].dx(0), u[0].dx(1), 0.], [u[1].dx(0), u[1].dx(1), 0.], [0., 0., u[1]/r]])

def div_cyl(u):
	return u[1]/r + u[0].dx(0) + u[1].dx(1)

# # NS equations Divergence form
# A1 =  eta*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx 
# #A1 =  eta*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*(inner(I,grad_cyl(w)))*r*dx 
# A2 = q*div_cyl(u)*r*dx 
# # Assembling of the system of eqs (4) + (5) + (6) in weak form
# A = A1 + A2


# Alternative form (from the tutorial to comprae with the iterative solution)
A1 = inner(grad_cyl(u), grad_cyl(w))*r*dx + div_cyl(w)*p*r*dx 
A2 = q*div_cyl(u)*r*dx
# Assembling of the system of eqs 
A = A1 + A2

f = Constant((0.0, 0.0))
L = inner(f, w)*r*dx

#L = - p_right*dot(n, w)*r*ds(5) 
#----------------------------------------------------------------------


# Assemble system
# A, bb = assemble_system(a, L, bcs)
U = Function(CHS)
solve(A == L, U, bc_tot, solver_parameters={'linear_solver':'lu'})


# Get sub-functions
p, u = U.split()

# Save solution in VTK format
ufile_pvd = File("velocity.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << p


def strain_rate(u): 
	return 0.5*(grad_cyl(u) + grad_cyl(u).T)

def sigma(u,p): 
	return -1.*p*I + 2*eta*strain_rate(u)

T = TensorFunctionSpace(mesh, 'CG', 1, (3,3) )

s_rate = strain_rate(u)
s = sigma(u,p)


strain_rate = project(s,T)
stress = project(s,T)

File("stress.pvd") << stress
File("strain_rate.pvd") << strain_rate




















