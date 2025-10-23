import numpy as np
import dolfinx
import mpi4py
import mpi4py.MPI
import ufl
import basix
from dolfinx.mesh import create_rectangle, CellType
from dolfinx import io, fem
import gmsh
import os
from mpi4py import MPI

def create_mesh(input_parameters):
    # Geometrical parameters
    r_c     = 0.90e-04
    r_tumor=input_parameters["r_tumor"]
    # Dimension of the problem
    gdim = 2
    # Initialize GMSH
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("Capsule")
    #Creation of the geometry of the mesh
    Tissue = gmsh.model.occ.addCircle(0, 0, 0, r_c,angle1=0.,angle2=np.pi/2, tag=-1,zAxis=[])
    #This function create a circle arc between angle 1 and angle 2, zAxis at the end is to define the normal of the circle.
    Tumor  = gmsh.model.occ.addCircle(0,0, 0, r_tumor,angle1=0.,angle2=np.pi/2, tag=-1,zAxis=[])
    # since this function create a line, we need to create all the other lines and the surface by hand
    O    = gmsh.model.occ.addPoint(0, 0, 0)
    lbx  = gmsh.model.occ.addLine(1,3)
    lsx  = gmsh.model.occ.addLine(3,O)
    lby  = gmsh.model.occ.addLine(2, 4)
    lsy  = gmsh.model.occ.addLine(4, O)
    cl1 = gmsh.model.occ.addCurveLoop([Tissue, lby, Tumor, lbx])
    s1  = gmsh.model.occ.addPlaneSurface([1], cl1)
    cl2 = gmsh.model.occ.addCurveLoop([lsx, lsy, 2])
    s2  = gmsh.model.occ.addPlaneSurface([2], cl2)
    # Use Fragment to split the rectangle and the circle, ensuring all boundaries are preserved
    # fragment = gmsh.model.occ.fragment([(2, Tissue)], [(2, Tumor)])
    gmsh.model.occ.synchronize()
    # Export the geometry to .geo_unrolled format for inspection
    # Remove existing file if it exists
    if os.path.exists("Capsule.geo_unrolled"):
        os.remove("Capsule.geo_unrolled")
    # Export the geometry
    gmsh.write("Capsule.geo_unrolled")
    # Define Physical Groups for boundary conditions
    tissue_limit, tag_bot, tag_left, tag_tumor = 1, 2, 3, 4
    # Create Physical Groups
    gmsh.model.addPhysicalGroup(gdim-1, [1], tissue_limit)
    gmsh.model.setPhysicalName(gdim-1, tissue_limit, 'tissue_limit')
    gmsh.model.addPhysicalGroup(gdim-1, [3,4], tag_bot)
    gmsh.model.setPhysicalName(gdim-1, tag_bot, 'bot_line')
    gmsh.model.addPhysicalGroup(gdim-1, [5,6], tag_left)
    gmsh.model.setPhysicalName(gdim-1, tag_left, 'left_line')
    gmsh.model.addPhysicalGroup(gdim-1, [2], tag_tumor)
    gmsh.model.setPhysicalName(gdim-1, tag_tumor, 'Tumor_line')
    gmsh.model.addPhysicalGroup(gdim, [1,2], 1)
    gmsh.model.setPhysicalName(gdim, 1, 'Tissue')
    # Set mesh parameters
    lc = 2.e-06
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
    # Generate and save the mesh
    gmsh.model.mesh.generate(gdim)
    # Ensure the GMSH model is properly created before writing
    gmsh.model.occ.synchronize()
    # Remove existing file if it exists
    if os.path.exists("Capsule.msh"):
        os.remove("Capsule.msh")
    gmsh.write("Capsule.msh")
    gmsh.finalize()
    #
    mesh, cell_markers, facet_tag = io.gmshio.read_from_msh('Capsule.msh', MPI.COMM_WORLD, 0, gdim)
    # Vérification des tags physiques détectés
    print("Tags disponibles dans facet_markers:", np.unique(facet_tag.values))
    #
    x = ufl.SpatialCoordinate(mesh)
    r = ufl.SpatialCoordinate(mesh)[0]
    normal   = ufl.FacetNormal(mesh)
    # Redefiniction dx and ds
    metadata = {"quadrature_degree": 4}
    dx = ufl.Measure("dx", domain = mesh, metadata=metadata)
    ds = ufl.Measure("ds", domain = mesh, subdomain_data = facet_tag)
    return mesh, facet_tag, normal, metadata, r, dx, ds
#
def create_capsule_mesh(input_parameters,output_dir):
    # caracteristic length
    lc = input_parameters[0]
    # Geometrical parameters
    r_int_caps     = input_parameters[1]
    r_ext_caps     = input_parameters[2]
    r_ext_tiss     = input_parameters[3]
    r_tumor        = input_parameters[4]
    #
    # Dimension of the problem
    gdim = 2
    #
    # Initialize GMSH
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("Confined_Capsule")
    #
    #This function create a circle arc between angle 1 and angle 2, zAxis at the end is to define the normal of the circle.
    Turmor      = gmsh.model.occ.addCircle(0,0, 0, r_tumor,angle1=0.,angle2=np.pi/2, tag=-1,zAxis=[])
    Int_Capsule = gmsh.model.occ.addCircle(0,0, 0, r_int_caps,angle1=0.,angle2=np.pi/2, tag=-1,zAxis=[])
    #Creation of the geometry of the mesh
    Ext_Capsule = gmsh.model.occ.addCircle(0, 0, 0, r_ext_caps,angle1=0.,angle2=np.pi/2, tag=-1,zAxis=[])
    Ext_Tissue  = gmsh.model.occ.addCircle(0, 0, 0, r_ext_tiss,angle1=0.,angle2=np.pi/2, tag=-1,zAxis=[])
    #
    # since this function create a line, we need to create all the other lines and the surface by hand
    #
    O               = gmsh.model.occ.addPoint(0, 0, 0)
    #
    ltumorx         = gmsh.model.occ.addLine(O,1)
    ltumor_capsulex = gmsh.model.occ.addLine(1,3)
    lcapsux         = gmsh.model.occ.addLine(3,5)
    ltissux         = gmsh.model.occ.addLine(5,7)
    #
    ltumory         = gmsh.model.occ.addLine(O,2)
    ltumor_capsuley = gmsh.model.occ.addLine(2,4)
    lcapsuy         = gmsh.model.occ.addLine(4,6)
    ltissuy         = gmsh.model.occ.addLine(6,8)
    #
    cl1 = gmsh.model.occ.addCurveLoop([ltumorx, Turmor, ltumory])
    s1  = gmsh.model.occ.addPlaneSurface([1], cl1)
    #
    cl2 = gmsh.model.occ.addCurveLoop([ltumor_capsulex, Int_Capsule, ltumor_capsuley,Turmor])
    s2  = gmsh.model.occ.addPlaneSurface([2], cl2)
    #
    cl3 = gmsh.model.occ.addCurveLoop([lcapsux, Ext_Capsule,lcapsuy, Int_Capsule])
    s3  = gmsh.model.occ.addPlaneSurface([3], cl3)
    #
    cl4 = gmsh.model.occ.addCurveLoop([ltissux, Ext_Tissue,ltissuy, Ext_Capsule])
    s4  = gmsh.model.occ.addPlaneSurface([4], cl4)
    gmsh.model.occ.synchronize()
    #
    # Export the geometry to .geo_unrolled format for inspection
    # Remove existing file if it exists
    if os.path.exists(f"{output_dir}/Confined_Capsule.geo_unrolled"):
        os.remove(f"{output_dir}/Confined_Capsule.geo_unrolled")
    #
    # Export the geometry
    gmsh.write(f"{output_dir}/Confined_Capsule.geo_unrolled")
    #
    #
    # Define Physical Groups for boundary conditions
    tag_symmetryx, tag_symmetryy,tag_tumor, tag_inner_shell, tag_Outer_shell, tag_Tissue_limit= 1, 2, 3, 4, 5, 6
    #
    # Create Physical Groups
    gmsh.model.addPhysicalGroup(gdim-1, [5,6,7,8], tag_symmetryx)
    gmsh.model.setPhysicalName(gdim-1, tag_symmetryx, 'symmetryx')
    gmsh.model.addPhysicalGroup(gdim-1, [9,10,11,12], tag_symmetryy)
    gmsh.model.setPhysicalName(gdim-1, tag_symmetryy, 'symmetryy')
    #
    gmsh.model.addPhysicalGroup(gdim-1, [1], tag_tumor)
    gmsh.model.setPhysicalName(gdim-1, tag_tumor, 'tumor_line')
    #
    gmsh.model.addPhysicalGroup(gdim-1, [2], tag_inner_shell)
    gmsh.model.setPhysicalName(gdim-1, tag_inner_shell, 'inner_shell')
    #
    gmsh.model.addPhysicalGroup(gdim-1, [3], tag_Outer_shell)
    gmsh.model.setPhysicalName(gdim-1, tag_Outer_shell, 'Outer_shell')
    #
    gmsh.model.addPhysicalGroup(gdim-1, [4], tag_Tissue_limit)
    gmsh.model.setPhysicalName(gdim-1, tag_Tissue_limit, 'Tissue_limit')
    #
    gmsh.model.addPhysicalGroup(gdim, [1,2,4], 1)
    gmsh.model.setPhysicalName(gdim, 1, 'Tissue')
    #
    gmsh.model.addPhysicalGroup(gdim, [3], 2)
    gmsh.model.setPhysicalName(gdim, 2, 'Capsule')
    #
    # gmsh.model.addPhysicalGroup(gdim, [2], 3)
    # gmsh.model.setPhysicalName(gdim, 3, 'Pericapsular')
    #
    # gmsh.model.addPhysicalGroup(gdim, [4], 4)
    # gmsh.model.setPhysicalName(gdim, 4, 'Tissue')
    # Set mesh parameters
    #
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
    #
    # Generate and save the mesh
    gmsh.model.mesh.generate(gdim)
    # Ensure the GMSH model is properly created before writing
    gmsh.model.occ.synchronize()
    #
    # Remove existing file if it exists
    if os.path.exists(f"{output_dir}/Confined_Capsule.msh"):
        os.remove(f"{output_dir}/Confined_Capsule.msh")
    gmsh.write(f"{output_dir}/Confined_Capsule.msh")
    gmsh.finalize()
    #
    # --- Read the mesh into dolfinx ---
    mesh, cell_markers, facet_markers = io.gmshio.read_from_msh(f"{output_dir}/Confined_Capsule.msh", MPI.COMM_WORLD, 0, gdim)
    #
    # Print detected boundary tags
    print("Available facet tags:", np.unique(facet_markers.values))
    #
    # --- Coordinates and integration measures ---
    x = ufl.SpatialCoordinate(mesh)
    r = x[0]  # radial coordinate
    #
    # Integration settings
    metadata = {"quadrature_degree": 10}
    dx = ufl.Measure("dx", domain=mesh, metadata=metadata)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_markers)
    normal = ufl.FacetNormal(mesh)
    return mesh, cell_markers,facet_markers,gdim, normal, metadata, r, dx, ds, lc, r_int_caps, r_tumor
def create_function_spaces(mesh):
    P1   = basix.ufl.element("Lagrange", "triangle", degree=1)
    CHDP = dolfinx.fem.functionspace(mesh, basix.ufl.mixed_element([P1, P1, P1, P1]))
    H    = dolfinx.fem.functionspace(mesh, P1)
    return CHDP, H

def readme_MB():
    return None