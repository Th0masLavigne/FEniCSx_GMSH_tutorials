import gmsh
from mpi4py import MPI
import numpy as np

# Initialize Gmsh
gmsh.initialize()
gmsh.option.setNumber('Geometry.Tolerance', 1e-4) 
mesh_comm = MPI.COMM_WORLD
model_rank = 0

#Characteristic length
lc = 5e-4  # The default larger mesh size away from the pipes
err= 1e-4
#===================================================================================
#          USEFULL DIMENSION OF THE DIFFERENT PARTS OF THE CHIP
#===================================================================================
#Chip body
[Length_b, Width_b, Height_b] = [35e-3,13e-3, 11e-3]
[Length, Width, Height] = [32e-3,1.5e-3, 6e-3]
#Wells
x1, y1, z1, dx1, dy1, dz1, r1 = 0, 0, 0, 0, 0, 3e-3, 7e-3#Main well where the exchange between air and medium will take place
x2, y2, z2, dx2, dy2, dz2, r2 = 0, 0, 0, 0, 0, -1e-3, 3e-3#Intermediary well that allow a better renewal of the culture medium
x3, y3, z3, dx3, dy3, dz3, r3 = 0, 0, -1e-3, 0, 0, -1e-3, 1.5e-3#Culture well where the cells will stay still
#Pipes
x4, y4, z4, dx4, dy4, dz4, r4 = -(21.5e-3), 0, 0.85e-3, 2*(21.5e-3), 0, 0, 0.85e-3#dimensions of the pipes, adapted for a 1/16" tubing
#Observation chambers
x5, y5, z5, dx5, dy5, dz5, r5 = 0, 0, dz1, 0, 0, 0.5e-3, 12.5e-3 #Holes dedicated to light and microscope obeservation
#Cover glasses
x6, y6, z6, dx6, dy6, dz6, r6 = 0, 0, 0, 0, 0, 0.17e-3, 11e-3 #Dimension of the cover glasses that will be used
gdim = 3
# Entry/ Exit of air & medium
x7, y7, z7, dx7, dy7, dz7, r7= 0, 0, 0, 0, 0, 0 ,2.2e-3 #Dimension of the 

axis_of_rotation = [dx4, dy4, dz4] 

def create_mesh(mesh, cell_type, prune_z=False):
    import meshio
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read": [cell_data]})
    return out_mesh
    
    
#===================================================================================
#                         CREATION OF THE MESH OUT OF VOLUMES
#===================================================================================
# Main code block where the geometry is set up and boundaries are categorized
if mesh_comm.rank == model_rank:
    #===================================================================================
    #                               CREATION OF AIR VOLUME
    #===================================================================================
    """========================CREATION OF THE BODY OF THE CHIP=============================="""
    # box2       = gmsh.model.occ.addBox(-(16e-3), 0, Height/2, Length , Width+1e-6, Height/2+2e-3,tag=-1)
    # """==============================CREATION OF THE WELLS==================================="""
    # AIR_Chamber  = gmsh.model.occ.addCylinder(x1, y1, 3e-3, dx1, dy1, 5e-3, r1,tag=-1,angle=np.pi)
    # New_chamber,_= gmsh.model.occ.cut([(3,AIR_Chamber)],[(3,box2)],removeTool=False,removeObject=True)
    """===============================CREATION OF THE PIPES=================================="""   
    E_air    = gmsh.model.occ.addCylinder(-15e-3, 0, 6e-3, 0, 0, 6.5e-3, r4, tag=-1, angle=np.pi)
    S_air     = gmsh.model.occ.addCylinder(15e-3, 0, 6e-3, 0, 0, 6.5e-3, r4, tag=-1, angle=np.pi)
    """========================= REDUCTION OF THE NUMBER OF VOLUMES=========================="""
    # Fragment_air,_= gmsh.model.occ.fragment(New_chamber,[(3,box2)])#[(3,E_air),(3,S_air),(3,box2)])
    # Total_air,_= gmsh.model.occ.fuse([Fragment_air[0]],[Fragment_air[x] for x in range(1,len(Fragment_air))])
    
    #===================================================================================
    #                         CREATION OF THE MEDIUM VOLUME
    #===================================================================================
    """========================CREATION OF THE BODY OF THE CHIP=============================="""
    box1       = gmsh.model.occ.addBox(-(16e-3), 0, 0, Length , Width, Height,tag=-1)
    """==============================CREATION OF THE WELLS==================================="""
    MED_Chamber = gmsh.model.occ.addCylinder(x1, y1, z1, dx1, dy1, Height-z1, r1,tag=-1,angle=np.pi)
    Well = gmsh.model.occ.addCylinder(x3, y3, z3, dx3, dy3, dz3, r3,tag=-1, angle=np.pi)
    IWell = gmsh.model.occ.addCylinder(x2, y2, z2, dx2, dy2, dz2, r2,tag=-1, angle=np.pi)
    """===============================CREATION OF THE PIPES=================================="""   
    Pipe_Left = gmsh.model.occ.addCylinder(-(21.5e-3), y4, 0.85e-3, 5.5e-3, dy4, dz4, r4, tag=-1, angle=np.pi)
    Pipe_Right = gmsh.model.occ.addCylinder((21.5e-3), y4, 0.85e-3, -5.5e-3, dy4, dz4, r4, tag=-1, angle=np.pi)
    gmsh.model.occ.rotate([(3, Pipe_Left),(3, Pipe_Right)], x4, y4, z4, *axis_of_rotation, np.pi)
    Fragment_fluid,_  = gmsh.model.occ.fragment([(3,MED_Chamber)],[(3,IWell),(3,Well),(3,box1),(3,E_air),(3,S_air),(3,Pipe_Left),(3,Pipe_Right)])
    """========================= REDUCTION OF THE NUMBER OF VOLUMES=========================="""
    Total_fluid,_= gmsh.model.occ.fuse([Fragment_fluid[0]],Fragment_fluid[1:])
    # Total_fragment,_= gmsh.model.occ.fragment(Total_fluid,Total_air)


    box1           = gmsh.model.occ.addBox(-(5e-2), -5e-2, -5e-2, 2*5e-2 , 2*5e-2, 5e-2+Height/2,tag=-1)
    Fragment_fluid,_ = gmsh.model.occ.intersect(Total_fluid,[(3,box1)],removeTool=True,removeObject=False)

    #===================================================================================
    #                         CREATION OF THE PLASTIC VOLUME
    #===================================================================================
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    # 
    Body_u   = gmsh.model.occ.addBox(-Length_b/2, 0, -(2.5e-3), Length_b, Width_b, Height_b,tag=-1)
    surfaces, volumes = [gmsh.model.getEntities(d) for d in [ 2, 3]]
    # 
    # """===============================CREATION OF THE PIPES==================================""" 
    E_air_p  = gmsh.model.occ.addCylinder(-15e-3, 0, 8.5e-3, 0, 0, 4e-3, r7, tag=-1, angle=np.pi)   
    S_air_p   = gmsh.model.occ.addCylinder( 15e-3, 0,8.5e-3, 0, 0, 4e-3, r7, tag=-1, angle=np.pi)
    E_medium_p = gmsh.model.occ.addCylinder(-(21.5e-3), 0, 8.5e-4, 4e-3, 0, 0, r7, tag=-1, angle=np.pi)
    S_medium_p = gmsh.model.occ.addCylinder( (21.5e-3), 0, 8.5e-4, -4e-3, 0,0 , r7, tag=-1, angle=np.pi)
    gmsh.model.occ.rotate([(3, E_medium_p),(3, S_medium_p)], x4, y4, z4, *axis_of_rotation, np.pi)
    plastic, _ = gmsh.model.occ.fuse([(3,Body_u)],[(3,E_air_p),(3,E_medium_p),(3,S_air_p),(3,S_medium_p)])
    # 
    """============================SET POSITION OF COVER GLASSES============================="""
    Obs_up = gmsh.model.occ.addCylinder(x5, y5, (8e-3-1e-3), dx5, dy5, dz5+1e-3, r5,tag=-1,angle=np.pi)
    
    Obs_down = gmsh.model.occ.addCylinder(x5, y5, -2.5e-3, dx5, dy5, dz5, r5,tag=-1,angle=np.pi)
    test2,_=gmsh.model.occ.cut(plastic, [(3, Obs_down)], removeTool=True, removeObject=False)    
    
    verre,_=gmsh.model.occ.cut(plastic,[(3,Obs_up)], removeTool=True, removeObject=False)
    # 
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    gmsh.write('Geo_oxy_3D.geo_unrolled')


    tdim =3 
    gmsh.model.addPhysicalGroup(tdim, [3], 1)
    gmsh.model.setPhysicalName(tdim, 1, 'air')
    gmsh.model.addPhysicalGroup(tdim, [2], 2)
    gmsh.model.setPhysicalName(tdim, 2, 'fluid')
    gmsh.model.addPhysicalGroup(tdim, [4], 3)
    gmsh.model.setPhysicalName(tdim, 3, 'plastic')
    gmsh.model.addPhysicalGroup(tdim, [5,6], 4)
    gmsh.model.setPhysicalName(tdim, 4, 'glass')


    # tdim =2 
    # gmsh.model.addPhysicalGroup(tdim, bottom, 3)
    # gmsh.model.setPhysicalName(tdim, 3, 'Bottom')
    # gmsh.model.addPhysicalGroup(tdim, symmetry, 4)
    # gmsh.model.setPhysicalName(tdim, 4, 'Symmetry')


    gmsh.option.setNumber("General.Terminal",1)
    gmsh.option.setNumber("Mesh.Optimize", True)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 5e-5)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 5e-4)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
    print("Min mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMin"))
    print("Max mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMax"))


    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write("Mesh.msh")
    gmsh.finalize()
