# Other GMSH Examples

This directory contains other examples of mesh generated with GMSH:
- a sphere from a stl file,
- a breast-like phantom as well as a tire-like geometry from *[Lavigne et al. 2023b](https://doi.org/10.1016/j.cma.2023.115889)*

The GMSH tagged meshes can be converted to meshes with fenics legacy compatibility considering the following:

```python
def create_mesh(mesh, cell_type, prune_z=False):
    import meshio
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={
                           "name_to_read": [cell_data]})
    return out_mesh
```
## 3D meshes:
```python
import meshio
msh = meshio.read("mesh.msh")
# BCs
surface_mesh = create_mesh(msh, "triangle", prune_z=False)
meshio.write("facet_tags.xdmf", surface_mesh)
# Mesh
tetra_mesh = create_mesh(msh, "tetra", prune_z=False)
meshio.write("mesh.xdmf", tetra_mesh)
```

## 2D meshes
```python
import meshio
msh = meshio.read("Mesh_refine.msh")
# BCs
line_mesh = create_mesh(msh, "line", prune_z=True)
meshio.write("facet_mesh_refine.xdmf", line_mesh)
# Mesh
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
meshio.write("mesh_refine.xdmf", triangle_mesh)
```

# Tomography to conform mutlipart meshes from several STL files
- *[Tomo2FE github](https://github.com/ANR-MultiFIRE/TomoToFE/blob/main/workflow2/Workflow2-Python.ipynb)*
- *[Tomo2FE article](https://letters.rilem.net/index.php/rilem/article/view/184)*

# Python tutorials
- *[Gitlab GMSH](https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/tutorials/python)*
- *[Remesh STL](https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t13.py)*
Contents (from *[Gitlab GMSH](https://gitlab.onelab.info/gmsh/gmsh/-/tree/master/tutorials/python)*):
* t1: Geometry basics, elementary entities, physical groups
* t2: Transformations, extruded geometries, volumes
* t3: Extruded meshes, parameters, options
* t4: Built-in functions, holes in surfaces, annotations, entity colors
* t5: Mesh sizes, loops, holes in volumes
* t6: Transfinite meshes, deleting entities
* t7: Background meshes
* t8: Post-processing, image export and animations
* t9: Plugins
* t10: Mesh size fields
* t11: Unstructured quadrangular meshes
* t12: Cross-patch meshing with compounds
* t13: Remeshing an STL file without an underlying CAD model
* t14: Homology and cohomology computation
* t15: Embedded points, lines and surfaces
* t16: Constructive Solid Geometry, OpenCASCADE geometry kernel
* t17: Anisotropic background mesh
* t18: Periodic meshes
* t19: Thrusections, fillets, pipes, mesh size from curvature
* t20: STEP import and manipulation, geometry partitioning
* t21: Mesh partitioning

Extended tutorials (API only):

* x1: Geometry and mesh data
* x2: Mesh import, discrete entities, hybrid models, terrain meshing
* x3: Post-processing data import: list-based
* x4: Post-processing data import: model-based
* x5: Additional geometrical data: parametrizations, normals, curvatures
* x6: Additional mesh data: integration points, Jacobians and basis functions
* x7: Additional mesh data: internal edges and faces