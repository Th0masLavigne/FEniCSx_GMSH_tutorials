# Thomas Lavigne
# 02-08-2024
# 
# From https://doi.org/10.1016/j.jmbbm.2023.105902
# 
#----------------------------------------------------------------------
# Libraries
#----------------------------------------------------------------------
# 
import dolfinx
import numpy
from dolfinx.io   import XDMFFile
from mpi4py       import MPI
#
# 
#----------------------------------------------------------------------
# Geometry
#----------------------------------------------------------------------
# 
## Box 
# Dimensions of the sample
[Length, Width, Height] = [6e-4, 2.5e-4, 4e-5]
# Discretization
[nx,ny,nz] = [30,15,8]
mesh = dolfinx.mesh.create_box(MPI.COMM_WORLD,numpy.array([[0.0,0.0,0.0],[Length, Width, Height]]), [nx,ny,nz], cell_type=dolfinx.mesh.CellType.tetrahedron)
def test_on_boundary(x):
	return numpy.logical_and((numpy.sqrt(numpy.power(x[0]-3e-4,2)+numpy.power(x[1],2))>=0.7*1.5e-4), (numpy.sqrt(numpy.power(x[0]-3e-4,2)+numpy.power(x[1],2))<=1.3*1.5e-4))
#
#
#----------------------------------------------------------------------
# Local Refinement
#----------------------------------------------------------------------
# 
refine_boudaries = [(11, lambda x: test_on_boundary(x))]
#
for _ in numpy.arange(2):
	# Refinement 
	refine_indices, refine_markers = [], []
	fdim = mesh.topology.dim-2
	for (marker, locator) in refine_boudaries:
		facets = dolfinx.mesh.locate_entities(mesh, fdim, locator)
		refine_indices.append(facets)
		refine_markers.append(numpy.full_like(facets, marker))
	refine_indices = numpy.hstack(refine_indices).astype(numpy.int32)
	refine_markers = numpy.hstack(refine_markers).astype(numpy.int32)
	# indices in meshtag must be sorted
	sorted_facets_refine = numpy.argsort(refine_indices)
	refine_tag = dolfinx.mesh.meshtags(mesh, fdim, refine_indices[sorted_facets_refine], refine_markers[sorted_facets_refine])
	mesh.topology.create_entities(fdim)
	mesh, _, _ = dolfinx.mesh.refine(mesh, refine_indices[sorted_facets_refine])
#
# 
#----------------------------------------------------------------------
# Mark the boundaries and domains
#----------------------------------------------------------------------
# 
def Omega_top(x):
	return numpy.logical_and((x[2] == Height), (numpy.sqrt(numpy.power(x[0]-3e-4,2)+numpy.power(x[1],2))>=1.39e-4))
#
def Omega_loading(x):
	return numpy.logical_and((x[2] == Height), (numpy.sqrt(numpy.power(x[0]-3e-4,2)+numpy.power(x[1],2))<=1.45e-4))
#
# It is complicated to have a clean circular boundary ... I recommend sticking to the GMSH local refinement applied before computations
# Create the facet tags (identify the boundaries)
# 1 = loading, 2 = top minus loading, 3 = bottom, 4 = left, 5 = right, 6 = Front, 7 = back
boundaries = [(1, lambda x: Omega_loading(x)),
              (2, lambda x: Omega_top(x)),
              (3, lambda x: numpy.isclose(x[2], 0.0)),
              (4, lambda x: numpy.isclose(x[0], 0.0)),
              (5, lambda x: numpy.isclose(x[0], Length)),
              (6, lambda x: numpy.isclose(x[1], 0.0)),
              (7, lambda x: numpy.isclose(x[1], Width))]
# Mark them
facet_indices, facet_markers = [], []
fdim = mesh.topology.dim - 1
for (marker, locator) in boundaries:
	facets = dolfinx.mesh.locate_entities_boundary(mesh, fdim, locator)
	facet_indices.append(facets)
	facet_markers.append(numpy.full_like(facets, marker))
facet_indices = numpy.hstack(facet_indices).astype(numpy.int32)
facet_markers = numpy.hstack(facet_markers).astype(numpy.int32)
sorted_facets = numpy.argsort(facet_indices)
facet_tag = dolfinx.mesh.meshtags(mesh, fdim, facet_indices[sorted_facets], facet_markers[sorted_facets])
facet_tag.name = "facets"
# Write XDMF
mesh.topology.create_connectivity(mesh.topology.dim-1, mesh.topology.dim)
with XDMFFile(mesh.comm, "facet_tags.xdmf", "w") as xdmftag:
    xdmftag.write_mesh(mesh)
    xdmftag.write_meshtags(facet_tag,mesh.geometry)
xdmftag.close()