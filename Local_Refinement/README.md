# Local Refinement

Local refinement can be perormed both in FEniCSx or in GMSH. Be aware that in case of complex geometries, a refinement based on a surface definition is more stable than a refinement from a square. This is for instance the case here with a cylinder inside of a cube (please refer to: *[Lavigne et al., 2023](https://doi.org/10.1016/j.jmbbm.2023.105902)* appendices where the method is described). In case of a choice of a local refinement inside of FEniCSx, **be careful with the facet tag definition** which might contain mismatches or not cover the whole surface (for instance with the boundary of the cylinder).

## Using FEniCSx

The idea for local refinement within FEniCSx is to identify the boundaries to be refined using a (marker,locator) tuple. 

```python
refine_boudaries = [(11, lambda x: test_on_boundary(x))]
```

Then, a loop corresponding to the number of operations we require is created to refine the mesh as many times as required:

```python
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
    mesh = dolfinx.mesh.refine(mesh, refine_indices[sorted_facets_refine])
```

## Using GMSH

The mesh refinement in GMSH relies on a distance field to an existing surface of the mesh. A precise tutorial is available [here](https://jsdokken.com/src/tutorial_gmsh.html).

```python3
# Create a distance field from which a local resolution is defined for the refinement
indenter_interface = surfaces[0][1]
distance = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance, "FacesList", [indenter_interface])
# A threshold function is defined:
resolution = r/10
threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
gmsh.model.mesh.field.setNumber(threshold, "LcMin", resolution)
gmsh.model.mesh.field.setNumber(threshold, "LcMax", 5*resolution)
gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.6*r)
gmsh.model.mesh.field.setNumber(threshold, "DistMax", r)
# If several fields are defined:
minimum = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [threshold]) # add other fields in the list if needed
gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
```