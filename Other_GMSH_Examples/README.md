# Other GMSH Examples

This directory contains other examples of mesh generated with GMSH:
- a sphere from a stl file,
- a breast-like phantom as well as a tire-like geometry from *[Lavigne et al. 2023b](https://doi.org/10.1016/j.cma.2023.115889)*

The GMSH tagged meshes can be converted to meshes with fenics legacy compatibility considering the following:
3D meshes:
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

2D meshes
```python
import meshio
msh = meshio.read("Mesh_refine.msh")
# BCs
line_mesh = create_mesh(msh, "line", prune_z=True)
meshio.write("facet_mesh_refine.xdmf", line_mesh)
# Mesh
triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
meshio.write("mesh_refine.xdmf", triangle_mesh)
return line_mesh, triangle_mesh
```