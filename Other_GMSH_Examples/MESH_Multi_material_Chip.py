import gmsh
from mpi4py import MPI
import numpy as np

# Initialize Gmsh
gmsh.initialize()
mesh_comm = MPI.COMM_WORLD
model_rank = 0

#Characteristic length
lc = 5e-4

#Air/medium chambers
[Length, Width, Height] = [9e-3, 5e-3, 7e-3]

#Observation Chambers
x1, y1, z1, dx1, dy1, dz1, r1 = 0, 0, 0, 0, 0, 3e-3, 4e-3
x2, y2, z2, dx2, dy2, dz2, r2 = 0, 0, 0, 0, 0, -1e-3, 3e-3
x3, y3, z3, dx3, dy3, dz3, r3 = 0, 0, -1e-3, 0, 0, -1e-3, 1.5e-3
x4, y4, z4, dx4, dy4, dz4, r4 = -(10e-3 + r1+Length+5e-3), 0, 0.8e-3, 2*(10e-3+r1+Length+5e-3), 0, 0, 0.8e-3


gdim = 3

def create_mesh(mesh, cell_type, prune_z=False):
    import meshio
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={
                           "name_to_read": [cell_data]})
    return out_mesh

# Main code block where the geometry is set up and boundaries are categorized
if mesh_comm.rank == model_rank:
    # Create the cylinders
    Chamber = gmsh.model.occ.addCylinder(x1, y1, z1, dx1, dy1, dz1, r1,tag=1,angle=np.pi)
    IChamber = gmsh.model.occ.addCylinder(x2, y2, z2, dx2, dy2, dz2, r2,tag=2, angle=np.pi)
    Well = gmsh.model.occ.addCylinder(x3, y3, z3, dx3, dy3, dz3, r3,tag=3, angle=np.pi)

    # Create Pipe
    Pipe = gmsh.model.occ.addCylinder(x4, y4, z4, dx4, dy4, dz4, r4, tag=4, angle=np.pi)
    gmsh.model.occ.rotate([(3,4)], x4, y4, z4, dx4, dy4, dz4, np.pi)
    
    gmsh.model.occ.cut([(3, Pipe)], [(3, Chamber)], removeTool=False)
    gmsh.model.occ.synchronize()
    
    # Create Air/medium chambers
    box1      = gmsh.model.occ.addBox(-(r1+10e-3+9e-3), 0, 0, Length, Width, Height)
    box2      = gmsh.model.occ.addBox((r1+10e-3), 0, 0, Length, Width, Height)
    gmsh.model.occ.cut([(3, Pipe)], [(3, box1)], removeTool=False)
    gmsh.model.occ.synchronize()
    
    gmsh.model.occ.cut([(3, Pipe)], [(3, box2)], removeTool=False)
    gmsh.model.occ.synchronize()

    gmsh.model.occ.synchronize()

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    
    
    # Retrieve the volumes and add them to physical groups
    surfaces, volumes = [gmsh.model.getEntities(dim) for dim in [2, 3]]

    # Define the physical groups for volumes
    for i, vol in enumerate(volumes):
        vol_id = vol[1]  # Get the volume ID
        gmsh.model.addPhysicalGroup(3, [vol_id])
        gmsh.model.setPhysicalName(3, vol_id, f'cyl{i+1}')

    #====================================================================================
    #                       THE PART UPSIDE THIS TEXT IS OK
    #====================================================================================
    # Categorize each boundary surface
    bot_marker, top_marker, Chamber_marker, IChamber_marker, Well_marker, bot_Chamber_marker,bot_IChamber_marker,Pipe_marker,Entry_marker = 1,2,3,4,5,6,7,8,9
    bot, top, mChamber, mIChamber, mWell, botChamber, botIChamber,Entry= [79],[68],[67,71,72],[73,76,77],[78,80,81],[69],[74],[49]
    Pipe=[38,40,42,43,45,47,48,50,52,57,58,60]
    
    rb_front_marker, rb_back_marker, rb_left_marker, rb_right_marker, rb_up_marker, rb_down_marker=10, 11, 12, 13, 14, 15
    lb_front_marker, lb_back_marker, lb_left_marker, lb_right_marker, lb_up_marker, lb_down_marker=16, 17, 18, 19, 20, 21
    rb_front, rb_back, rb_left, rb_right, rb_up, rb_down=[84],[85],[82],[83],[87],[86]
    lb_front, lb_back, lb_left, lb_right, lb_up, lb_down=[66,54,56],[64],[61],[62],[63],[65]
    boundaries = gmsh.model.getBoundary(volumes, oriented=False)
    # for boundary in boundaries:
        # recup.append(boundary[1])

    # Assign the categorized surfaces to physical groups and set their names
    # ADD CULTURE AND OBSERVATION CHAMBERS
    gmsh.model.addPhysicalGroup(boundaries[0][0], bot, bot_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], bot_marker, 'bot')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], top, top_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], top_marker, 'top')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], mChamber, Chamber_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], Chamber_marker, 'mChamber')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], mIChamber, IChamber_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], IChamber_marker, 'mIChamber')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], mWell, Well_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], Well_marker, 'mWell')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], botChamber, bot_Chamber_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], bot_Chamber_marker, 'botChamber')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], botIChamber, bot_IChamber_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], bot_IChamber_marker, 'botIChamber')
    
    #ADD PIPES
    gmsh.model.addPhysicalGroup(boundaries[0][0], Pipe, Pipe_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], Pipe_marker, 'Pipe')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], Entry, Entry_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], Entry_marker, 'Entry')
    
    #ADD AIRE/MEDIUM CHAMBERS
    gmsh.model.addPhysicalGroup(boundaries[0][0], rb_front, rb_front_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], rb_front_marker, 'rb_front')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], rb_back, rb_back_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], rb_back_marker, 'rb_back')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], rb_left, rb_left_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], rb_left_marker, 'rb_left')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], rb_right, rb_right_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], rb_right_marker, 'rb_right')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], rb_up, rb_up_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], rb_up_marker, 'rb_up')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], rb_down, rb_down_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], rb_down_marker, 'rb_down')
    
    
    
    
    
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], lb_front, lb_front_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], lb_front_marker, 'lb_front')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], lb_back, lb_back_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], lb_back_marker, 'lb_back')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], lb_left, lb_left_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], lb_left_marker, 'lb_left')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], lb_right, lb_right_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], lb_right_marker, 'lb_right')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], lb_up, lb_up_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], lb_up_marker, 'lb_up')
    
    gmsh.model.addPhysicalGroup(boundaries[0][0], lb_down, lb_down_marker)
    gmsh.model.setPhysicalName(boundaries[0][0], lb_down_marker, 'lb_down')
    
   
    # gmsh.model.addPhysicalGroup(boundaries[0][0], recup, recup_marker)
    # gmsh.model.setPhysicalName(boundaries[0][0], recup_marker, 'recup')
    
    gmsh.model.occ.synchronize()
    gmsh.write('Geo_FLUIDIQUE.geo_unrolled')

    gmsh.model.occ.synchronize()
    
    #===================================================================================
    #                           REFINMENT OF THE MESH
    #===================================================================================
    pipe_surfaces = []
    for s in surfaces:
        center = gmsh.model.occ.getCenterOfMass(s[0], s[1])
        # Assuming the Pipe's surfaces have distinct z-coordinates
        # Replace with appropriate condition to identify Pipe's surfaces
        if center[2] == z4: 
            pipe_surfaces.append(s[1])
    
    # Begin adaptive mesh refinement setup
    distance_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance_field, "SurfacesList", pipe_surfaces)

    resolution_min = 1e-3  # The desired smaller mesh size on the pipes
    resolution_max = lc  # The default larger mesh size away from the pipes
    threshold_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold_field, "IField", distance_field)
    gmsh.model.mesh.field.setNumber(threshold_field, "LcMin", resolution_min)
    gmsh.model.mesh.field.setNumber(threshold_field, "LcMax", resolution_max)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 2*r4)  # Start decreasing mesh size at this distance
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", 5*r4)  # End decreasing mesh size at this distance

    # Apply the threshold field as the mesh size field
    minimum_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum_field, "FieldsList", [threshold_field])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum_field)

    # Synchronize the model before meshing
    gmsh.model.occ.synchronize()

    # Now, generate the mesh
    gmsh.model.mesh.generate(3)
    
    gmsh.write("mesh.msh")
    gmsh.finalize()

    
import meshio
msh = meshio.read("mesh.msh")
# BCs
surface_mesh = create_mesh(msh, "triangle", prune_z=False)
meshio.write("facet_mesh_refined.xdmf", surface_mesh)
# Mesh
tetra_mesh = create_mesh(msh, "tetra", prune_z=False)
meshio.write("mesh.xdmf", tetra_mesh)
    
exit()