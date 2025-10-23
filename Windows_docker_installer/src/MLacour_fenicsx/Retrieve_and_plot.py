#Retrieve the values of saturation, capillay pressure for poromechanical model
import dolfinx
import basix.ufl
import ufl
import random as random
import numpy as np
import os
import gmsh
import sys
import petsc4py
import matplotlib.pyplot as plt
#
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, io, plot
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.fem.petsc import NonlinearProblem
# --- Extract DOF coordinates and pressure values ---
"""
collapsed_aggregate_pressure = _pressure.collapse()
dof_coords_pressure = collapsed_aggregate_pressure.function_space.tabulate_dof_coordinates()
pressure_values = collapsed_aggregate_pressure.x.array

if dof_coords_pressure.ndim == 1:
    dof_coords_pressure = dof_coords_pressure.reshape((-1, collapsed_aggregate_pressure.function_space.mesh.geometry.dim))


# --- Calculate and extract Cell Saturation values ---
# Get the cell fraction component from YN
cell_fraction_ufl = YN.sub(1)

# Create a Function to interpolate the UFL expression for cell saturation
V_scalar = W2 # Assuming W2 is a suitable scalar FunctionSpace



cell_saturation_func = fem.Function(V_scalar)

epsilon_sat = (cell_fraction_ufl - epsilon_s)
filtered_sat = 0.5 * abs(epsilon_sat) + 0.5 * epsilon_sat
cell_saturation_func.interpolate(fem.Expression(filtered_sat / (1 - epsilon_s), V_scalar.element.interpolation_points()))

# Direct access: No need to collapse, as it's already a standalone Function
dof_coords_cell_saturation = cell_saturation_func.function_space.tabulate_dof_coordinates()
cell_saturation_values = cell_saturation_func.x.array

if dof_coords_cell_saturation.ndim == 1:
    dof_coords_cell_saturation = dof_coords_cell_saturation.reshape((-1, cell_saturation_func.function_space.mesh.geometry.dim))

# --- Calculate and extract Capillary Pressure values ---
# Create a Function to interpolate the UFL expression for capillary pressure
capillary_pressure_func = fem.Function(V_scalar)
# Corrected: Use cell_saturation_func instead of cell_saturation_0 and V_scalar instead of W1
capillary_pressure_func.interpolate(fem.Expression(ufl.conditional(ufl.lt(r, r_int_caps), a_func * ufl.tan(ufl.pi * cell_saturation_func / 2), 0.0), V_scalar.element.interpolation_points()))

# Direct access: No need to collapse, as it's already a standalone Function
dof_coords_capillary_pressure = capillary_pressure_func.function_space.tabulate_dof_coordinates()
capillary_pressure_values = capillary_pressure_func.x.array

# New addition: Ensure values outside the capsule are explicitly zero for plotting
if dof_coords_capillary_pressure.ndim == 1:
    dof_coords_capillary_pressure = dof_coords_capillary_pressure.reshape((-1, capillary_pressure_func.function_space.mesh.geometry.dim))

# Calculate radial distance for each DOF point for capillary pressure
distances_capillary = np.linalg.norm(dof_coords_capillary_pressure[:, :2] - center, axis=1)

# Create a mask for points that are geometrically outside the active region.
# Using r_int_caps - 1e-11 consistent with your conditional.
outside_capsule_mask = distances_capillary > (r_int_caps)

# This is the line that performs the action you suggested:
capillary_pressure_values[outside_capsule_mask] = 0.0


# --- Capsule quarter-circle arc ---
theta = np.linspace(0, np.pi / 2, 200)
x_arc = r_int_caps * np.cos(theta)
y_arc = r_int_caps * np.sin(theta)

# --- Plotting with Subplots ---
fig, axes = plt.subplots(1, 3, figsize=(24, 8)) # 1 row, 3 columns for 3 plots

# Plot 1: Aggregate Pressure
ax0 = axes[0]
sc0 = ax0.scatter(dof_coords_pressure[:, 0], dof_coords_pressure[:, 1],
                  c=pressure_values, cmap='coolwarm', s=15)
plt.colorbar(sc0, ax=ax0, label="Aggregate Pressure (Pa)")
ax0.plot(x_arc, y_arc, '--', linewidth=2, color='blue', label='Capsule Boundary (¼-circle)')
ax0.set_aspect("equal")
ax0.set_xlabel("X (m)")
ax0.set_ylabel("Y (m)")
ax0.set_title("Aggregate Pressure Field")
ax0.legend()
ax0.grid(True)

# Plot 2: Cell Saturation
ax1 = axes[1]
sc1 = ax1.scatter(dof_coords_cell_saturation[:, 0], dof_coords_cell_saturation[:, 1],
                  c=cell_saturation_values, cmap='viridis', s=15)
plt.colorbar(sc1, ax=ax1, label="Cell Saturation")
ax1.plot(x_arc, y_arc, '--', linewidth=2, color='blue', label='Capsule Boundary (¼-circle)')
ax1.set_aspect("equal")
ax1.set_xlabel("X (m)")
ax1.set_ylabel("Y (m)")
ax1.set_title("Cell Saturation Field")
ax1.legend()
ax1.grid(True)

# Plot 3: Capillary Pressure
ax2 = axes[2]
sc2 = ax2.scatter(dof_coords_capillary_pressure[:, 0], dof_coords_capillary_pressure[:, 1],
                  c=capillary_pressure_values, cmap='plasma', s=15)
plt.colorbar(sc2, ax=ax2, label="Capillary Pressure (Pa)")
ax2.plot(x_arc, y_arc, '--', linewidth=2, color='blue', label='Capsule Boundary (¼-circle)')
ax2.set_aspect("equal")
ax2.set_xlabel("X (m)")
ax2.set_ylabel("Y (m)")
ax2.set_title("Capillary Pressure Field")
ax2.legend()
ax2.grid(True)

# Adjust layout and save the figure
fig.tight_layout()
fig.savefig(f"{output_dir}Combined_Pressure_Saturation_Capillary.png")
plt.close(fig)

print("✅ Saved: 'Combined_Pressure_Saturation_Capillary.png'")

"""

"""
==========================================================================================================

                                        PLOTTING FUNCTIONS

==========================================================================================================
"""
def evaluate_point(mesh, function, contributing_cells, point, output_list, index, variable_name="Variable"):
    """
    Evaluates a function at a given point and stores the result in an output list.

    Parameters:
        mesh: Computational mesh.
        function: The function to be evaluated.
        contributing_cells: List of mesh cells contributing to the evaluation.
        point: The point where the function is evaluated.
        output_list: The list to store the output value.
        index: The index at which to store the value.
        variable_name (str, optional): Name of the variable being evaluated.
    """
    if contributing_cells.size > 0:
        function_eval = function.eval(point, contributing_cells[:1])
        if function_eval is None:
            return
        #
        if mesh.comm.rank == 0:
            result = mesh.comm.gather(function_eval, root=0)
            value = next((res[0] for res in result if res is not None), 0.0)
            output_list[index] = value

# Evolution of values along the X-axis

def plot_x_axis_evolution(x_values, data_sets, ylabel, filename,output_dir):
    # Plots the evolution of a variable along the X-axis at three specific time steps
    # as subplots within a single figure.

    # Args:
        # x_values (np.array): The x-coordinates (distance).
        # data_sets (list of tuples): A list of (index, data_array, label_string) for each time snapshot.
                                     # Expected to be in order: (t0, t_mid, t_end).
        # ylabel (str): Label for the Y-axis.
        # filename (str): Name of the file to save the plot (e.g., "Displacement_Evolution_vs_Distance.jpg").

    # Create a figure with 3 subplots arranged vertically, sharing the X-axis for easy comparison
    fig, axes = plt.subplots(3, 1, figsize=(10, 15), sharex=True)
    fig.suptitle(f"{ylabel} along X-axis Evolution", fontsize=16) # Main title for the entire figure

    # Iterate through data_sets (which should be t0, tmid, tend) and plot each on its own subplot
    for i, (index_val, data, label) in enumerate(data_sets):
        ax = axes[i] # Get the current subplot axis (ax[0] for t0, ax[1] for tmid, ax[2] for tend)
        ax.plot(x_values, data, label=f"Time: {label}", color='blue') # Consistent color for now
        ax.set_ylabel(ylabel) # Each subplot gets its own Y-axis label
        ax.set_title(label) # Use the time label as the subplot title
        ax.grid(True, linestyle="--", linewidth=0.5)
        ax.legend(loc='best') # Add legend to each subplot

    axes[-1].set_xlabel("Distance along X-axis (m)") # Only the bottom subplot needs the X-axis label

    fig.tight_layout(rect=[0, 0.03, 1, 0.96]) # Adjust layout to prevent suptitle overlap with subplots
    fig.savefig(os.path.join(output_dir, filename))
    plt.close(fig)


    
def plot_combined_x_axis_evolutions(plot_data_list, main_title, output_filename, output_dir):
    """
    Plots multiple variables' evolution along the X-axis at three specific time steps
    as subplots within a single figure (2x2 grid).

    Args:
        plot_data_list (list): A list of dictionaries, where each dictionary
                               contains data for one variable to be plotted:
                               - 'x_values': The x-coordinates (distance).
                               - 'data_sets': A list of (index, data_array, label_string)
                                              for each time snapshot (t0, t_mid, t_end).
                               - 'ylabel': Label for the Y-axis of this variable's subplot.
                               - 'title': Title for this variable's subplot.
        main_title (str): Main title for the entire figure.
        output_filename (str): Name of the file to save the plot.
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), sharex=True)
    axes = axes.flatten()  # Flatten the 2x2 array of axes for easy iteration

    fig.suptitle(main_title, fontsize=18)

    for i, plot_info in enumerate(plot_data_list):
        ax = axes[i]
        x_values = plot_info['x_values']
        data_sets = plot_info['data_sets']
        ylabel = plot_info['ylabel']
        title = plot_info['title']

        for _, data, label in data_sets:
            ax.plot(x_values, data, label=f"Time: {label}")

        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle="--", linewidth=0.5)
        ax.legend(loc='best')

    # Set common X-axis label for the bottom row
    for i in [2, 3]: # Indices for the bottom row subplots
        axes[i].set_xlabel("Distance along X-axis (m)")

    fig.tight_layout(rect=[0, 0.03, 1, 0.96]) # Adjust layout to prevent suptitle overlap
    fig.savefig(os.path.join(output_dir, output_filename))
    plt.close(fig)
    print(f"✅ Combined plot saved: {os.path.join(output_dir, output_filename)}")

    
def plot_all_center_evolutions(str_list,time_values, pressure_center, saturation_center, oxygen_center, potential_center, output_dir):
    fig, axs = plt.subplots(4, 1, figsize=(8, 12), sharex=True)
    title, legend1,  legend2, legend3, legend4 = str_list[0], str_list[1], str_list[2], str_list[3], str_list[4],
    # Plot Pressure
    axs[0].plot(np.array(time_values) / 3600, pressure_center[:len(time_values)], color='blue')
    axs[0].set_ylabel(legend1)
    axs[0].legend([legend1])
    axs[0].grid(True, linestyle="--", linewidth=0.5)

    # Plot Saturation
    axs[1].plot(np.array(time_values) / 3600, saturation_center[:len(time_values)], color='orange')
    axs[1].set_ylabel(legend2)
    axs[1].legend([legend2])
    axs[1].grid(True, linestyle="--", linewidth=0.5)

    # Plot Oxygen
    axs[2].plot(np.array(time_values) / 3600, oxygen_center[:len(time_values)], color='green')
    axs[2].set_ylabel(legend3)
    axs[2].legend([legend3])
    axs[2].grid(True, linestyle="--", linewidth=0.5)

    # Plot Potential
    axs[3].plot(np.array(time_values) / 3600, potential_center[:len(time_values)], color='red')
    axs[3].set_ylabel(legend4)
    axs[3].set_xlabel("Time (h)")
    axs[3].legend([legend4])
    axs[3].grid(True, linestyle="--", linewidth=0.5)

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{title}Evolution_Center_All.jpg"))
    plt.close(fig)
    
    
# Define the new setup function
def setup_simulation_data_and_analysis_points(mesh, lc, tf, dt, output_dir, target_point=None):
    """
    Sets up simulation data arrays, finds target points/cells, and performs initial mesh visualization.

    Args:
        mesh: The Dolfinx mesh object.
        lc (float): Characteristic length, used for tolerance.
        tf (float): Total simulation time.
        dt (float): Initial time step.
        output_dir (str): Directory for output files.
        target_point (np.array, optional): The desired coordinate for the center point.
                                          Defaults to np.array([15e-5, 15e-5, 0.]).

    Returns:
        tuple: A tuple containing:
            - time_values (list): List to store time points.
            - heaviside_center (list): List to store Heaviside values at the center.
            - point_center (np.array): The closest mesh point to the target_point.
            - colliding_cells_center (dolfinx.geometry.BBTree): Cells colliding with point_center.
    """
    if target_point is None:
        target_point = np.array([15e-5, 15e-5, 0.])

    # Extract mesh coordinates
    all_points = mesh.geometry.x[:]

    # Find points near the X-axis (y ≈ 0)
    tolerance = lc / 2
    points_x = all_points[np.abs(all_points[:, 1]) < tolerance]
    points_x = points_x[np.argsort(points_x[:, 0])]
    num_points_x = len(points_x)
    print(f"✅ Number of points detected on the X-axis: {num_points_x}")

    # Detect corresponding mesh cells for points_x (this part might not be used directly later for 'center' point)
    tree = dolfinx.geometry.bb_tree(mesh, mesh.geometry.dim)
    cell_candidates_x = dolfinx.geometry.compute_collisions_points(tree, points_x)
    colliding_cells_x = dolfinx.geometry.compute_colliding_cells(mesh, cell_candidates_x, points_x)

    # Define storage arrays for time evolution
    num_steps = int(tf / dt) + 1 # Note: This num_steps is an initial estimate. time_values and heaviside_center are lists to handle variable dt.
    time_values = []
    heaviside_center = []

    # Define indices for three key time steps (useful for reference, though num_steps is dynamic here)
    index_t0, index_tmid, index_tend = 0, num_steps - 1, num_steps - 1
    print("\n index_t0 ", index_t0, "index_tmid ", index_tmid, "index_tend ", index_tend)

    # Compute Euclidean distance and retrieve closest point to target_point
    distances = np.linalg.norm(all_points - target_point, axis=1)
    closest_idx = np.argmin(distances)
    point_center = all_points[closest_idx]
    print(f"✅ Closest selected point to {target_point}: {point_center}")

    # Detect corresponding mesh cells for the center point
    cell_candidates_center = dolfinx.geometry.compute_collisions_points(tree, point_center.reshape(1, -1))
    colliding_cells_center = dolfinx.geometry.compute_colliding_cells(mesh, cell_candidates_center, point_center.reshape(1, -1))

    # 🔹 Visualization of the mesh and the selected point
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(all_points[:, 0], all_points[:, 1], s=5, color='gray', label="Mesh Points")
    ax.scatter(target_point[0], target_point[1], s=100, color='red', marker='x', label="Selected Target Point")
    ax.scatter(point_center[0], point_center[1], s=100, color='blue', marker='o', label="Closest Retained Point")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_title("Visualization of Selected Point vs Mesh")
    ax.legend()
    ax.grid(True)
    plot_path = os.path.join(output_dir, "mesh_point_visualization.png")
    fig.savefig(plot_path)
    plt.close(fig)
    print("✅ Visualization image saved as 'mesh_point_visualization.png'")

    return time_values, heaviside_center, point_center, colliding_cells_center
#
# Define the new function for evaluating and storing Heaviside values
def evaluate_and_store_function(point_center, colliding_cells_center, H1, t, function_center_list):
    """
    Evaluates the H (H1) function at a specific center point,
    prints its value, and stores it in a list.

        point_center (np.array): The coordinates of the center point [x, y, z].
        colliding_cells_center (dolfinx.geometry.BBTree.compute_colliding_cells result):
                                 The result from computing colliding cells for the center point.
        H1 (dolfinx.fem.Function or ufl.Expression): The Heaviside function to evaluate.
        t (float): The current simulation time.
        heaviside_center_list (list): The list to which the evaluated H1_value will be appended.

    Returns:
        float or None: The evaluated Heaviside (H1) value at the current time step, or None if no cell is found.
    """
    H1_value = None
    cells_center = colliding_cells_center.links(0) # Retrieve contributing cells associated with point_center

    if len(cells_center) > 0:
        # Evaluate H1 (H function) at the target point
        H1_value = H1.eval(point_center.reshape(1, -1), np.array([cells_center[0]]))[0]
    else:
        # This branch correctly handles the case where no cell is found
        print(f"⚠️ No contributing cell found for point {point_center} at time {t}.")

    # Verify and display values
    print(f"📌 Time {t}: {H1}={H1_value}")

    # Store evaluated Heaviside value
    if H1_value is not None:
        function_center_list.append(H1_value)

    return H1_value

def plot_function_center_evolution(time_values, function_values, title, output_filename, output_dir, ylabel_text="Function Value"):
    """
    Plots the evolution of a function at a specific point over time.

    Args:
        time_values (list): List of time points.
        function_values (list): List of function values corresponding to time_values.
        title (str): Title of the plot.
        output_filename (str): Name of the file to save the plot.
        output_dir (str): Directory where the plot will be saved.
        ylabel_text (str, optional): Text for the Y-axis label. Defaults to "Function Value".
    """
    plt.figure(figsize=(10, 6))
    plt.plot(time_values, function_values, marker='o', linestyle='-', color='purple')
    plt.title(title)
    plt.xlabel("Time (s)")
    plt.ylabel(ylabel_text) # Use the new ylabel_text parameter
    plt.grid(True)
    plt.tight_layout()
    plot_path = os.path.join(output_dir, output_filename)
    plt.savefig(plot_path)
    plt.close()
    print(f"✅ Plot saved: '{plot_path}'")

def plot_lame_coeff():
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    # ===============================
    #   VISUALIZATION SECTION
    # ===============================
    # 1. Compute cell centroids for plotting (geometric centers of mesh cells)
    tdim = mesh.topology.dim                                # Topological dimension of the mesh (2D or 3D)
    mesh.topology.create_connectivity(tdim, 0)              # Ensure connectivity between cells and vertices is built
    cell_vertices = mesh.topology.connectivity(tdim, 0)     # Get cell-to-vertex connectivity
    x = mesh.geometry.x                                     # Get coordinates of all mesh vertices
    #
    # Compute centroids of each cell by averaging its vertex coordinates
    cell_centroids = np.array([np.mean(x[cell_vertices.links(i)], axis=0) for i in range(mesh.topology.index_map(tdim).size_local)])
    #
    # 2. Set up colormaps and normalization for lambda and mu values
    lambda_vals = np.unique(lambda_array)   # Unique values of lambda (used for discrete coloring)
    mu_vals = np.unique(mu_array)           # Unique values of mu
    #
    # Assign color to each unique value
    lambda_colors = ["purple", "gold", "blue", "red"][:len(lambda_vals)]
    mu_colors = ["navy", "orange", "green", "brown"][:len(mu_vals)]
    #
    # Create colormaps from selected colors
    lambda_cmap = ListedColormap(lambda_colors)
    mu_cmap = ListedColormap(mu_colors)
    #
    # Create normalization boundaries for discrete color mapping
    lambda_norm = BoundaryNorm(lambda_vals.tolist() + [lambda_vals[-1] + 1], len(lambda_vals))
    mu_norm = BoundaryNorm(mu_vals.tolist() + [mu_vals[-1] + 1], len(mu_vals))
    #
    # 3. Plotting the spatial distribution of lambda and mu
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)  # Create side-by-side subplots
    #
    # --- Lambda plot ---
    sc0 = axs[0].scatter(cell_centroids[:, 0], cell_centroids[:, 1], c=lambda_array, cmap=lambda_cmap, norm=lambda_norm, s=12)
    axs[0].set_title("Spatial Distribution of Lambda")
    axs[0].set_xlabel("X (m)")
    axs[0].set_ylabel("Y (m)")
    axs[0].set_aspect("equal")
    axs[0].grid(True)
    #
    # Add colorbar for lambda
    cbar0 = plt.colorbar(sc0, ax=axs[0], ticks=lambda_vals)
    cbar0.ax.set_yticklabels([f"{v:.0f}" for v in lambda_vals])  # Show values as integers
    cbar0.set_label("Lambda Value")
    #
    # --- Mu plot ---
    sc1 = axs[1].scatter(cell_centroids[:, 0], cell_centroids[:, 1], c=mu_array, cmap=mu_cmap, norm=mu_norm, s=12)
    axs[1].set_title("Spatial Distribution of Mu")
    axs[1].set_xlabel("X (m)")
    axs[1].set_ylabel("Y (m)")
    axs[1].set_aspect("equal")
    axs[1].grid(True)
    #
    # Add colorbar for mu
    cbar1 = plt.colorbar(sc1, ax=axs[1], ticks=mu_vals)
    cbar1.ax.set_yticklabels([f"{v:.0f}" for v in mu_vals])
    cbar1.set_label("Mu Value")
    #
    # Save figure to file
    plt.savefig("lambda_mu_discrete_colormap.png")
    plt.close()
    print("✅ Saved: 'lambda_mu_discrete_colormap.png'")
    return

def readme_RaP():
    return None
