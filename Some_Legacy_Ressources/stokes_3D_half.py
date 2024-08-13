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
#Making a cylindrical geometry
cylinder1 = Cylinder(Point(0, 0, 0), Point(0, 0, L1), H1, H1)
cylinder2 = Cylinder(Point(0, 0, L1), Point(0, 0, (L1 + L2)), H2, H2)
cube1     = Box(Point(H1, -H1, 0), Point(-H1, 0, (L1 + L2)))
cube2     = Box(Point(0, 0, 0), Point(-H1, H1, (L1 + L2)))

geometry = cylinder1 + cylinder2 - cube1 - cube2

# Making Mesh (40 corresponds to the mesh density)
mesh = generate_mesh(geometry, 40)

# Save the mesh
file_mesh = File('cylinder.pvd')
file_mesh << mesh
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of functional space, trial & test functions
#----------------------------------------------------------------------
#___________________________________________________________________________
# Define function space
CP1  = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
CV2  = VectorElement("Lagrange", mesh.ufl_cell(), 3)
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
tol  = 3.e-2
tolc = 3.e-2

# Bounds needed for the fluid flow problem

# Boundary_1(Left bound)
class Boundary_1(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[2], 0, tol)

bound_1 = Boundary_1()
bound_1.mark(boundary_markers, 1)



# Boundary_2(right bound)
class Boundary_2(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[2], (L1 + L2), tol)

bound_2 = Boundary_2()
bound_2.mark(boundary_markers, 2)


# Boundary_3U(upper circular bound)
xc1 = 0. 
yc1 = 0.
class Boundary_3U(SubDomain):
	def inside(self, x, on_boundary):
		return abs(sqrt((x[1]-yc1)*(x[1]-yc1)+(x[0]-xc1)*(x[0]-xc1)) - H1) < tolc and on_boundary

bound_3U = Boundary_3U()
bound_3U.mark(boundary_markers, 3)


# Boundary_4U(upper circular bound small)
xc1 = 0. 
yc1 = 0.
class Boundary_4U(SubDomain):
	def inside(self, x, on_boundary):
		return abs(sqrt((x[1]-yc1)*(x[1]-yc1)+(x[0]-xc1)*(x[0]-xc1)) - H2) < tolc and x[2] >= L1 and on_boundary

bound_4U = Boundary_4U()
bound_4U.mark(boundary_markers, 4)

# Boundaries writting in file
file_boundary = File('boundary.pvd')
file_boundary << boundary_markers



# Boundary_5 (intermediate plan)
xc1 = 0. 
yc1 = 0.
class Boundary_5(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[2], L1, tol) and sqrt((x[1]-yc1)*(x[1]-yc1)+(x[0]-xc1)*(x[0]-xc1)) >= (H2-tol) and on_boundary

bound_5 = Boundary_5()
bound_5.mark(boundary_markers, 5)


# Boundary_6 (horizontal symmetry plan)
class Boundary_6(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0, tol)

bound_6 = Boundary_6()
bound_6.mark(boundary_markers, 6)


# Boundary_7 vertical symmetry plan)
class Boundary_7(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], 0, tol)

bound_7 = Boundary_7()
bound_7.mark(boundary_markers, 7)



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

v_null   = Constant((0.0, 0.0, 0.0))
v_in     = Constant((0.0, 0.0, 1.0))

zero     = Constant(0.)
p_right  = Constant(0.)

# No-slip boundary condition for velocity
bcNS3    = DirichletBC(CHS.sub(1), v_null, boundary_markers, 3)
bcNS4    = DirichletBC(CHS.sub(1), v_null, boundary_markers, 4)
bcNS5    = DirichletBC(CHS.sub(1), v_null, boundary_markers, 5)

# in fow vx = 1
bcNS1 = DirichletBC(CHS.sub(1), v_in, boundary_markers, 1)

# out flow vy = 0
bcNS2 = DirichletBC(CHS.sub(1).sub(1), zero, boundary_markers, 2)

# horizontal symmetry plan
bcNS6 = DirichletBC(CHS.sub(1).sub(1), zero, boundary_markers, 6)

# vertical symmetry plan
bcNS7 = DirichletBC(CHS.sub(1).sub(0), zero, boundary_markers, 7)

# Collect boundary conditions
bc_tot = [bcNS1, bcNS2, bcNS3, bcNS4, bcNS5, bcNS6, bcNS7]
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Weak form
#----------------------------------------------------------------------
# Define variational problem
import dolfin as df
I = df.Identity(3)

# NS equations Divergence form
A1 =  eta*inner(grad(u), grad(w))*dx +  eta*inner(grad(u).T, grad(w))*dx -  p*(inner(I,grad(w)))*dx 
A2 = q*div(u)*dx 

# Assembling of the system of eqs (4) + (5) + (6) in weak form
A = A1 + A2

L = - p_right*dot(n, w)*ds(2) 
#----------------------------------------------------------------------


# Assemble system
# A, bb = assemble_system(a, L, bcs)
U = Function(CHS)
solve(A == L, U, bc_tot, solver_parameters={'linear_solver':'lu'})

# solve(A == L, U, bc_tot,solver_parameters={'linear_solver': 'gmres','preconditioner': 'ilu'})

# prm = parameters.krylov_solver # short form
# prm.absolute_tolerance = 1E-10
# prm.relative_tolerance = 1E-6
# prm.maximum_iterations = 1000

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




















