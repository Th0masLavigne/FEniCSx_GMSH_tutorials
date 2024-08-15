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

## Tomography to conform mutlipart meshes from several STL files
- *[Tomo2FE github](https://github.com/ANR-MultiFIRE/TomoToFE/blob/main/workflow2/Workflow2-Python.ipynb)*
- *[Tomo2FE article](https://letters.rilem.net/index.php/rilem/article/view/184)*