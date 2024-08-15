# Thomas Lavigne
# 01/31/2024 
#
# Useful jorgen dokken links + forum
# https://fenicsproject.discourse.group/
# https://jsdokken.com/tutorials.html
#
# 
# 
# 
def extract_coordinates_lines_line_loops(file_path):
    """
    Read a .geo file and returns the ID + corrdinates of the points and the connectivity table (Lines)
    """
    coordinates = {}
    lines = []
    line_loops = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip().startswith('Point'):
                parts = line.split()
                point_id = int(parts[0].strip('Point()'))
                coordinates[point_id] = [float(parts[3].rstrip(',')), float(parts[4].rstrip(',')), float(parts[5].rstrip(','))]
            elif line.strip().startswith('Line'):
                parts = line.split()
                line_nodes = [int(parts[3].strip('{,')),int(parts[4].strip('};'))]
                lines.append(line_nodes)
            elif line.strip().startswith('Line Loop'):
                parts = line.split()
                line_loop_nodes = [int(node_id) for node_id in parts[2:]]
                line_loops.append(line_loop_nodes)
    return coordinates, lines, line_loops



def create_mesh(mesh, cell_type, prune_z=False):
    import meshio
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={
                           "name_to_read": [cell_data]})
    return out_mesh

def create_gmsh_from_geo(file_path, lc, lc_refine, d_min, d_max, cadre, tag_plane=1, tag_left=2, tag_right=3, tag_top=4, tag_bottom=5, tag_inclusions=6):
    """
    Finir par le rectangle lorsque l'on crée le .geo
    Les tags sont définis au début
    lc: longueur charactéristique au niveau des points
    En cas de besoin, on peut rafiner localement autour des inclusions
    """
    import gmsh
    import numpy as np
    # 
    #
    if cadre == 0: # means first
        condition = 0
    elif cadre == 1: # measn last
        condition = 1
    else:
        print(".geo not respecting the required order") 
        exit()
    # Some kernel parameters of GMSH
    gmsh.initialize()
    gmsh.clear()
    gmsh.option.setNumber("General.Terminal",1)
    gmsh.option.setNumber("Mesh.Optimize", True)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", True)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.01*lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.0)
    print("Min mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMin"))
    print("Max mesh size allowed ", gmsh.option.getNumber("Mesh.MeshSizeMax")) 
    # gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    # gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    # gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # 
    gmsh.model.occ.synchronize()
    #
    # Apply the above function
    coordinates, lines, __ = extract_coordinates_lines_line_loops(file_path)
    #
    # Create the model
    gmsh.model.add("Chip")
    #
    # Identification of the lineloops, i.e., surfaces, from the connectivity
    lineloops_indices = []
    index_begin = lines[0][0]
    N_loops = 0
    line0 = lines[0]
    #
    for idx in coordinates:
        # Create the geometrical nodes
        gmsh.model.occ.addPoint(coordinates[idx][0],  coordinates[idx][1],  lc)
    for index in range(len(lines)):
        # Create the geometrical lines
        gmsh.model.occ.addLine(lines[index][0], lines[index][1])
        gmsh.model.occ.synchronize()
        if lines[index][1]==line0[0]:
            index_end = lines[index][0]
            lineloops_indices.append([index_begin, index_end])
            N_loops+=1
            if index < len(lines)-1:
                index_begin = lines[index+1][0]
                line0=lines[index+1]
    # 
    # Create the surfaces
    surf={} 
    cl={}       
    for x in range(N_loops):
        cl[x] = gmsh.model.occ.addCurveLoop([i for i in range(lineloops_indices[x][0],lineloops_indices[x][1]+1)])
        surf[x] = gmsh.model.occ.addPlaneSurface([x+1], cl[x])
    #
    gmsh.model.occ.synchronize()
    # identify the geometrical elements
    edges, surfaces, volumes = [gmsh.model.getEntities(d) for d in [1, 2, 3]]
    # Write a geo_unrolled for debug
    gmsh.write('Geom_beforecut_refine.geo_unrolled')
    #
    # Remove the inclusions by boolean operation
    if condition == 1:
        for ii in range(len(surfaces)-1):
            gmsh.model.occ.cut([surfaces[len(surfaces)-1]], [surfaces[ii]], tag = -1, removeObject = True, removeTool = True)
    else:
        for ii in range(len(surfaces)-1):
            gmsh.model.occ.cut([surfaces[0]], [surfaces[ii+1]], tag = -1, removeObject = True, removeTool = True)
    #
    # Remove duplicate entities and synchronize
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    #
    # update the geometrical elements
    points, edges, surfaces, volumes = [gmsh.model.getEntities(d) for d in [0,1, 2, 3]]
    gmsh.model.occ.synchronize()
    # Surface tagging
    tdim = 2 
    su = []
    xmin,ymin,zmin,xmax,ymax,zmax=1e10,1e10,1e10,-1e10,-1e10,-1e10
    # if condition==1:
    for iik in range(len(surfaces)):
        su.append(surfaces[iik][1])
    gmsh.model.addPhysicalGroup(tdim, su, tag_plane)
    for elem in su:
        x,y,z,xx,yy,zz=gmsh.model.getBoundingBox(tdim, elem)
        xmin=min(xmin,x)
        ymin=min(ymin,y)
        zmin=min(zmin,z)
        xmax=max(xmax,xx)
        ymax=max(ymax,yy)
        zmax=max(zmax,zz)
    # xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(tdim, su[-1])
    # else:
    #     gmsh.model.addPhysicalGroup(tdim, [surfaces[0][1]], tag_plane)
    #     xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(tdim, surfaces[0][1])
    print(xmin, ymin, zmin, xmax, ymax, zmax)
    gmsh.model.setPhysicalName(tdim, tag_plane, 'Mesh')
    gmsh.model.occ.synchronize()
    # Edges tagging
    bottom, top, left, right, inclusions = [],[],[],[],[]
    boundaries = gmsh.model.getBoundary(surfaces, oriented=False)
    print(edges==boundaries)
    for boundary in boundaries:
        center_of_mass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
        # if np.isclose(center_of_mass[0], min([coordinates[x][0] for x in coordinates])):
        if np.isclose(center_of_mass[0], xmin, rtol=1e-05, atol=1e-06):
            left.append(boundary[1])
        elif np.isclose(center_of_mass[0], xmax, rtol=1e-05, atol=1e-06):
            right.append(boundary[1])
        # elif np.isclose(center_of_mass[1], min([coordinates[y][1] for y in coordinates])):
        elif np.isclose(center_of_mass[1], ymin, rtol=1e-05, atol=1e-06):
            bottom.append(boundary[1])
        elif np.isclose(center_of_mass[1], ymax, rtol=1e-05, atol=1e-06): 
            top.append(boundary[1])
        else:
            inclusions.append(boundary[1])
    # 
    tdim = 1 
    # we associate a tag with the corresponding list of edges
    gmsh.model.addPhysicalGroup(tdim, left, tag_left)
    gmsh.model.setPhysicalName(tdim, tag_left, 'left')
    gmsh.model.addPhysicalGroup(tdim, top, tag_top)
    gmsh.model.setPhysicalName(tdim, tag_top, 'top')
    gmsh.model.addPhysicalGroup(tdim, right, tag_right)
    gmsh.model.setPhysicalName(tdim, tag_right, 'right')
    gmsh.model.addPhysicalGroup(tdim, bottom, tag_bottom)
    gmsh.model.setPhysicalName(tdim, tag_bottom, 'bottom')
    gmsh.model.addPhysicalGroup(tdim, inclusions, tag_inclusions)
    gmsh.model.setPhysicalName(tdim, tag_inclusions, 'inclusions')
    gmsh.model.occ.synchronize()
    # 
    # Write a geo_unrolled for debug
    gmsh.write('Geom_refine.geo_unrolled')
    # 
    # Local refinement : https://jsdokken.com/src/tutorial_gmsh.html
    distance = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(distance, "EdgesList", inclusions)
    # 
    threshold = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold, "IField", distance)
    gmsh.model.mesh.field.setNumber(threshold, "LcMin", lc_refine)
    gmsh.model.mesh.field.setNumber(threshold, "LcMax", lc)
    gmsh.model.mesh.field.setNumber(threshold, "DistMin", d_min)
    gmsh.model.mesh.field.setNumber(threshold, "DistMax", d_max)
    # 
    # If several local refinement, replicate the above procedure and complete the list
    minimum = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(
        minimum, "FieldsList", [threshold])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    gmsh.model.occ.synchronize()
    # 
    # Generate the mesh in 2D
    gmsh.model.mesh.generate(2)
    gmsh.write("Mesh_refine.msh")
    gmsh.finalize()
    # 
    # Convert to xdmf
    import meshio
    msh = meshio.read("Mesh_refine.msh")
    # BCs
    line_mesh = create_mesh(msh, "line", prune_z=True)
    meshio.write("facet_mesh_refine.xdmf", line_mesh)
    # Mesh
    triangle_mesh = create_mesh(msh, "triangle", prune_z=True)
    meshio.write("mesh_refine.xdmf", triangle_mesh)
    return line_mesh, triangle_mesh
    
if __name__ == "__main__":
    file_path = "chip5x5.geo"  # Chemin vers votre fichier gmsh.geo
    # Define tags
    tag_plane=1
    tag_left=2
    tag_right=3
    tag_top=4
    tag_bottom=5
    tag_inclusions=6
    # 
    #
    # Characteristic length
    lc = 1e-4
    lc_refine = lc/10
    d_min = 5e-5
    d_max = 5e-4
    #
    # cadre 0 if first, 1 if last
    cadre = 1
    line_mesh, triangle_mesh=create_gmsh_from_geo(file_path, lc, lc_refine, d_min, d_max, cadre, tag_plane, tag_left, tag_right, tag_top, tag_bottom, tag_inclusions)







