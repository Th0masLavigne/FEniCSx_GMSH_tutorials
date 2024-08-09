from dolfinx import default_scalar_type
from dolfinx.fem import (Constant, dirichletbc, Function, functionspace, assemble_scalar,
                         form, locate_dofs_geometrical, locate_dofs_topological)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile, gmshio
from dolfinx.mesh import create_unit_square, locate_entities
from dolfinx.plot import vtk_mesh

from ufl import (SpatialCoordinate, TestFunction, TrialFunction,
                 dx, grad, inner)

from mpi4py import MPI

import meshio
import gmsh
import numpy as np
import pyvista

pyvista.start_xvfb()

mesh = create_unit_square(MPI.COMM_WORLD, 10, 10)
Q = functionspace(mesh, ("DG", 0))

def Omega_0(x):
    return x[1] <= 0.5


def Omega_1(x):
    return x[1] >= 0.5

kappa = Function(Q)
cells_0 = locate_entities(mesh, mesh.topology.dim, Omega_0)
cells_1 = locate_entities(mesh, mesh.topology.dim, Omega_1)

kappa.x.array[cells_0] = np.full_like(cells_0, 1, dtype=default_scalar_type)
kappa.x.array[cells_1] = np.full_like(cells_1, 0.1, dtype=default_scalar_type)

V = functionspace(mesh, ("Lagrange", 1))
u, v = TrialFunction(V), TestFunction(V)
a = inner(kappa * grad(u), grad(v)) * dx
x = SpatialCoordinate(mesh)
L = Constant(mesh, default_scalar_type(1)) * v * dx
dofs = locate_dofs_geometrical(V, lambda x: np.isclose(x[0], 0))
bcs = [dirichletbc(default_scalar_type(1), dofs, V)]

from dolfinx import log
log.set_log_level(log.LogLevel.INFO)

problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

# Filter out ghosted cells
tdim = mesh.topology.dim
num_cells_local = mesh.topology.index_map(tdim).size_local
marker = np.zeros(num_cells_local, dtype=np.int32)
cells_0 = cells_0[cells_0 < num_cells_local]
cells_1 = cells_1[cells_1 < num_cells_local]
marker[cells_0] = 1
marker[cells_1] = 2
mesh.topology.create_connectivity(tdim, tdim)
topology, cell_types, x = vtk_mesh(mesh, tdim, np.arange(num_cells_local, dtype=np.int32))

# p = pyvista.Plotter(window_size=[800, 800])
# grid = pyvista.UnstructuredGrid(topology, cell_types, x)
# grid.cell_data["Marker"] = marker
# grid.set_active_scalars("Marker")
# p.add_mesh(grid, show_edges=True)
# if pyvista.OFF_SCREEN:
#     figure = p.screenshot("subdomains_structured.png")
# 
# p.show()

# p2 = pyvista.Plotter(window_size=[800, 800])
# grid_uh = pyvista.UnstructuredGrid(*vtk_mesh(V))
# grid_uh.point_data["u"] = uh.x.array.real
# grid_uh.set_active_scalars("u")
# p2.add_mesh(grid_uh, show_edges=True)
# if not pyvista.OFF_SCREEN:
#     p2.show()
# else:
#     figure = p2.screenshot("subdomains_structured2.png")

def eval_kappa(x):
    values = np.zeros(x.shape[1], dtype=default_scalar_type)
    # Create a boolean array indicating which dofs (corresponding to cell centers)
    # that are in each domain
    top_coords = x[1] > 0.5
    bottom_coords = x[1] < 0.5
    values[top_coords] = np.full(sum(top_coords), 0.1)
    values[bottom_coords] = np.full(sum(bottom_coords), 1)
    return values

kappa2 = Function(Q)
kappa2.interpolate(eval_kappa)

# Difference in kappa's
error = mesh.comm.allreduce(assemble_scalar(form((kappa - kappa2)**2 * dx)))
print('errer',error)