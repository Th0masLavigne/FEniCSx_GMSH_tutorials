# FEniCSx & GMSH Tutorial Repository (v0.9.0)

This repository contains all materials, scripts, and documentation for a Finite Element Modeling workshop, focusing on **FEniCSx v0.9.0** and **GMSH v4.11+**. The project is directly related to Thomas Lavigne's AFR-FNR research, a collaborative effort involving the University of Luxembourg, Bordeaux University (I2M), and ENSAM Paris (IBHGC); and to the workshop conducted at I2M Bordeaux of the 9th-11th September 2024.

The primary goal is to guide users through advanced finite element methods using **FEniCSx** and robust, complex mesh generation with **GMSH**.

Please give a look here <https://th0maslavigne.github.io/FEniCSx_GMSH_tutorials/README.html>

---

## 1. Repository Structure Overview

The repository is logically partitioned into **Installation** methodologies and practical **Tutorials**.

| Directory | Content Summary |
| :--- | :--- |
| **`Install/`** | Scripts and comprehensive documentation for setting up the FEniCSx environment on various platforms: Windows (WSL2, Docker) and High-Performance Computing (HPC) systems (Spack). |
| **`Tutoriels/`** | Python scripts demonstrating a range of finite element simulations (Elasticity, Fluid Dynamics, Thermal Analysis) and advanced GMSH meshing techniques. |

### Core Installation Sub-Folders Breakdown

| Sub-Folder | Target Platform | Primary Files |
| :--- | :--- | :--- |
| `Install/WSL_Docker_installer` | Windows (WSL2 + Docker) | `install_wsl2_fenicsx_docker.bat`, `FEniCSx_interactive.bat` |
| `Install/WSL_Local_installer` | Windows (WSL2 Native) | `install_local_fenicsx.bat`, `apt-installed-list.txt` |
| `Install/Windows_docker_installer`| Windows (Docker Desktop Direct) | `FEniCSx_docker_0-9-0.bat`, `Interactive_FEniCSx_0-9-0.bat` |
| `Install/Spack_Install_FEniCSx` | Linux / HPC | `Spack_INSTALL.md`, `Slurm_spack_fenicsx_0_9_0.sh` |

---

## 2. Requirements & Environment Setup

To successfully execute the provided examples, the following components are necessary:
1.  **FEniCSx v0.9.0** runtime environment (Requires one of the installation methods outlined below).
2.  **GMSH** software (version **>4.11**) for mesh generation.
3.  **ParaView** for post-processing and visualization.

### Installation Guidance

**Detailed, step-by-step instructions for all setup methods are located in the `Install/Install_processes.md` file.**

For Windows users, the **WSL2 Docker Installer** (`Install/WSL_Docker_installer/install_wsl2_fenicsx_docker.bat`) is the **most recommended** approach, offering a self-contained, high-performance Linux environment.

Here-after are links to the official download pages:
- [FEniCSx v0.9.0](https://fenicsproject.org/download/)
- [Docker Desktop](https://desktop.docker.com/win/main/amd64/165256/Docker%20Desktop%20Installer.exe) v4.34.0 or later *(macOS or linked with WSL2 on Windows)* 
- [Docker Engine](https://docs.docker.com/engine/install/ubuntu/) *(macOS or Linux or Ubuntu WSL2)*,
- [GMSH](https://gmsh.info/#Download) software,
- [Paraview](https://www.paraview.org/download/) software.

---

## 3. Usage and Execution Launchers 

### Environment Aliases (WSL2 Docker Method)

After running the **WSL2 Docker Installer**, custom aliases are available inside the resulting WSL distribution (`FEniCSxenv`):

| Alias | Function | Description |
| :--- | :--- | :--- |
| `fenicsx` | **Interactive Shell** | Launches an interactive Bash shell inside the FEniCSx container. |
| `fenicsx-run` | **Run Script/Command** | Executes a single command or Python script non-interactively. |
| `fenicsx-jupyter` | **Jupyter Lab** | Launches Jupyter Lab, exposing it for access via your web browser. |

### Note on Graphical Interfaces (Crucial) 

By default, Docker containers lack a graphical environment. Any script that attempts to open a GUI window (e.g., GMSH's interactive viewer or certain Matplotlib backends) within the container may fail or crash.

To prevent this when running GMSH-related scripts, you **must pass the `close` argument** to suppress the GUI:

```bash
# Example of running a GMSH script non-interactively
fenicsx-run "python3 Tutoriels/Other_GMSH_Examples/chip5x5.geo close"
````

-----


## 4\. FEniCSx & GMSH Workshop Contents

The workshop is designed to provide a comprehensive introduction to meshing with **GMSH** and finite element analysis with **FEniCSx**.

---

### 4.1. Finite Element Method (FEM) Overview

* A brief reminder about the **Finite Element Method** fundamentals.

---

### 4.2. GMSH Mesh Creation Techniques

The workshop explores versatile methods for generating complex meshes, focusing on control and export compatibility:

* **Geometry Definition:**
    * Creation from a **sketch**.
    * Creation from **elementary geometries**.
    * Utilizing **symmetries** and **boolean operations**.
    * Generating meshes from a pre-existing **STL file**.
* **Advanced Meshing:**
    * Implementing **local refinement procedures** (referencing Lavigne et al., 2023).
    * **Boundary and domain tagging** for material and boundary condition mapping.
* **Export:**
    * Export compatibility: **FEniCSx** vs. **Fenics Legacy**.

---

### 4.3. FEniCSx Computation and Problem Solving

The computational portion of the workshop covers a range of physics and solution methods:

* **Problem Types:**
    * **Thermal Problems:** Stationary and transient analysis.
    * **Fluid Dynamics:** **Stokes Equation** solving in **2D and 3D**.
    * **Poro-Mechanics:** **Terzaghi poromechanical model** (referencing Lavigne et al., 2023).
    * **Solid Continuous Mechanics:**
        * Elastic and **hyper-elastic** problems.
        * Problems involving **penalty contact**.
        * **Incremental** resolution (for updated Lagrangian formulation).
        * Evaluation of a quantity (e.g., reaction force).
* **Numerical Methods:**
    * **Linear and Non-Linear resolutions**.
    * **Updated mesh resolution** (mesh update). 

>**For updated lagrangian please contact Carla Cornillon**

---

### 4.4. Key Academic References and Resources

The specialized topics (local refinement, poromechanics, and incremental resolution) are supported by the author's published work, available in full on their GitHub repositories:

* **Porous Media and Local Refinement:**
    * **Lavigne et al., 2023** (DOI: [10.1016/j.jmbbm.2023.105902](https://doi.org/10.1016/j.jmbbm.2023.105902))
    * This work is related to the **Terzaghi poromechanical model** and **local refinement** in the workshop.
    * **GitHub Repository:** [Dolfinx Porous Media](https://github.com/Th0masLavigne/Dolfinx_Porous_Media.git)

* **Incremental Resolution and Mesh Update:**
    * **Lavigne et al., 2024** (Work presented for an example of **incremental resolution** and **mesh update**.)
    * This includes the use of volume tags for material mapping and **reaction force evaluation**.
    * **GitHub Repository:** [Skin Porous Modelling](https://github.com/Th0masLavigne/Skin_porous_modelling.git)

-----

## 5\. Utility Code Snippet

A simple Python function is provided here for exporting data lists to a CSV file:

```python
def export_to_csv(data, filename, header=None):
    import csv
    try:
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            if header:
                writer.writerow(header)
            writer.writerows(data)
        print(f"Data exported to {filename} successfully")
    except Exception as e:
        print(f"An error occurred while exporting data to {filename}: {e}")
```

-----

## 6\. Acknowledgments and Licensing

### Research and Funding

This activity is a core component of **Thomas Lavigne's PhD work**. This research was entirely or partially funded by the **Luxembourg National Research Fund (FNR)**, under grant reference No. 17013182.

### Open Access License

For the purpose of open access, the author has applied a **Creative Commons Attribution 4.0 International (CC BY 4.0)** license to this work.

### Inspiration and Support

The author extends thanks to:

  * **Jørgen S. Dokken** and **Christophe Geuzaine**, whose exemplary work and tutorials served as inspiration for this repository.
  * **F. Daghia** and **Ludovic Chamoin** for the high-quality finite element courses delivered at the École Normale Supérieure Paris-Saclay.
  * **Jack Hale** and **Stéphane Urcun** for their invaluable assistance in code debugging throughout the PhD journey.
  * **Giuseppe Sciumè** for the invitation to co-organize the workshop at the Bordeaux University.


-----

## 7\. Comprehensive Resource Links

| Category | Resource Title | Description / Context | Link |
| :--- | :--- | :--- | :--- |
| **FEniCSx** | FEniCS Project Website | Official project homepage and documentation portal. | [https://fenicsproject.org/](https://fenicsproject.org/) |
| **FEniCSx** | FEniCSx Tutorial | Comprehensive tutorial by Jørgen S. Dokken. | [https://jsdokken.com/dolfinx-tutorial/](https://jsdokken.com/dolfinx-tutorial/) |
| **FEniCSx** | FEniCSx Workshop (Deep Dive) | Workshop materials by Jørgen S. Dokken. | [https://jsdokken.com/FEniCS-workshop/src/deep\_dive/expressions.html](https://jsdokken.com/FEniCS-workshop/src/deep_dive/expressions.html) |
| **FEniCSx** | FEniCSx Discourse Forum | Official community forum for questions and support. | [https://fenicsproject.discourse.group/](https://fenicsproject.discourse.group/) |
| **FEniCSx** | FEniCSx GitHub Repositories | Main GitHub organization for the FEniCS project. | [https://github.com/orgs/FEniCS/repositories](https://github.com/orgs/FEniCS/repositories) |
| **FEniCSx** | DOLFINx Documentation (v0.8.0) | Core documentation for the Python interface. | [https://docs.fenicsproject.org/dolfinx/v0.8.0/python/](https://docs.fenicsproject.org/dolfinx/v0.8.0/python/) |
| **FEniCSx** | BASIX Documentation (v0.8.0) | Documentation for the finite element basis function library. | [https://docs.fenicsproject.org/basix/v0.8.0/python/](https://docs.fenicsproject.org/basix/v0.8.0/python/) |
| **FEniCSx** | UFL Documentation | Unified Form Language documentation. | [https://fenics.readthedocs.io/projects/ufl/en/latest/](https://fenics.readthedocs.io/projects/ufl/en/latest/) |
| **FEniCSx** | FFCx Documentation | FEniCS Form Compiler documentation. | [https://docs.fenicsproject.org/ffcx/main/](https://docs.fenicsproject.org/ffcx/main/) |
| **FEniCSx** | Solver Options (PETSc) | Manual for the PETSc library used for linear and nonlinear solvers. | [https://petsc.org/main/manual/ksp/](https://petsc.org/main/manual/ksp/) |
| **FEniCSx** | Changelog | Release notes for DOLFINx on GitHub. | [https://github.com/FEniCS/dolfinx/releases](https://github.com/FEniCS/dolfinx/releases) |
| **FEniCSx** | FEniCS Legacy Documentation | Documentation for older versions of FEniCS. | [https://fenicsproject.org/olddocs/](https://fenicsproject.org/olddocs/) |
| **FEniCSx** | FEniCS Book | The classic FEniCS book (2011). | [https://launchpadlibrarian.net/83776282/fenics-book-2011-10-27-final.pdf](https://launchpadlibrarian.net/83776282/fenics-book-2011-10-27-final.pdf) |
| **FEniCSx** | Mesh Update Discussion | Discourse thread on incrementally deformed mesh schemes/mesh geometry. | [https://fenicsproject.discourse.group/t/how-to-do-updated-lagrangian-when-the-displacement-lives-in-a-different-space-to-the-mesh-geometry/10760/2](https://fenicsproject.discourse.group/t/how-to-do-updated-lagrangian-when-the-displacement-lives-in-a-different-space-to-the-mesh-geometry/10760/2) |
| **FEniCSx** | Parallel Run Discussion | Discourse thread on discrepancies between serial and parallel runs. | [https://fenicsproject.discourse.group/t/different-results-in-serial-and-parallel-run-dolfinx/4370](https://fenicsproject.discourse.group/t/different-results-in-serial-and-parallel-run-dolfinx/4370) |
| **FEniCSx** | Parallel Demo (PML) | Example of a parallel run with Perfectly Matched Layers. | [https://github.com/FEniCS/dolfinx/blob/main/python/demo/demo\_pml.py](https://github.com/FEniCS/dolfinx/blob/main/python/demo/demo_pml.py) |
| **FEniCSx** | Parallel Demo (Axisymmetry) | Example of a parallel run for an axisymmetric problem. | [https://github.com/FEniCS/dolfinx/blob/main/python/demo/demo\_axis.py](https://github.com/FEniCS/dolfinx/blob/main/python/demo/demo_axis.py) |
| **GMSH** | GMSH Download Page | Official download page for GMSH software. | [https://gmsh.info/\#Download](https://www.google.com/search?q=https://gmsh.info/%23Download) |
| **GMSH** | GMSH Manual | Detailed documentation and reference manual. | [https://gmsh.info/doc/texinfo/gmsh.html](https://gmsh.info/doc/texinfo/gmsh.html) |
| **GMSH** | GMSH Introductory Presentation | General overview course material (PDF). | [https://gmsh.info/doc/course/general\_overview.pdf](https://gmsh.info/doc/course/general_overview.pdf) |
| **GMSH** | GMSH GitLab Repository | Source code and issue tracker. | [https://gitlab.onelab.info/gmsh/gmsh](https://gitlab.onelab.info/gmsh/gmsh) |
| **GMSH** | GMSH API Tutorials | Tutorials on using the GMSH API. | [https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api/](https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/api/) |
| **GMSH** | OpenCascade Commands | Reference for OpenCascade commands within GMSH. | [https://koehlerson.github.io/gmsh.jl/dev/occ/occ/](https://koehlerson.github.io/gmsh.jl/dev/occ/occ/) |
| **Post-processing**| ParaView Download Page | Official download page for the visualization tool. | [https://www.paraview.org/download/](https://www.paraview.org/download/) |
| **Docker** | Docker Website | Official website for Docker Desktop and Docker Engine. | [https://www.docker.com/products/docker-desktop/](https://www.docker.com/products/docker-desktop/) |
| **Docker** | Docker Desktop Installer (Win) | Direct link to the Windows Docker Desktop installer. | [Docker Desktop Installer](https://www.google.com/search?q=https://desktop.docker.com/win/main/amd64/165256/Docker%2520Desktop%2520Installer.exe) |
| **Docker** | Docker Desktop Installer (Mac ARM64) | Direct link for macOS with Apple Silicon. | [Docker for Mac (ARM64)](https://www.google.com/search?q=https://desktop.docker.com/mac/main/arm64/165256/Docker.dmg) |
| **Docker** | Docker Desktop Installer (Mac AMD64) | Direct link for macOS with Intel chips. | [Docker for Mac (AMD64)](https://www.google.com/search?q=https://desktop.docker.com/mac/main/amd64/165256/Docker.dmg) |
| **Docker** | Docker Installation Resource | Official guide for Windows installation. | [https://docs.docker.com/desktop/install/windows-install/](https://docs.docker.com/desktop/install/windows-install/) |
| **Docker** | Docker Hub | Repository for Docker images. | [https://hub.docker.com/](https://hub.docker.com/) |
| **Docker** | Docker Cheat Sheet (Official) | Official quick reference guide. | [https://docs.docker.com/get-started/docker\_cheatsheet.pdf](https://docs.docker.com/get-started/docker_cheatsheet.pdf) |
| **Docker** | Docker Cheat Sheet (Collabnix) | Alternative detailed cheat sheet. | [https://dockerlabs.collabnix.com/docker/cheatsheet/](https://dockerlabs.collabnix.com/docker/cheatsheet/) |
| **Docker** | Windows Docker Cheat Sheet | Gist for Windows-specific Docker commands. | [https://gist.github.com/danijeljw/a7a2553bd06742648172363ce3983a9a](https://gist.github.com/danijeljw/a7a2553bd06742648172363ce3983a9a) |
| **Containers** | Singularity (Apptainer) | Documentation for the HPC-focused container platform. | [https://docs.sylabs.io/guides/3.5/user-guide/introduction.html](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) |
| **HPC** | Spack Package Manager | Official GitHub repository for the Spack package manager. | [https://github.com/spack/spack/](https://github.com/spack/spack/) |
| **HPC Containers** | Spack Containers Intro | Quick introduction to using Spack for creating container images (Docker/Singularity). | [https://spack.readthedocs.io/en/latest/containers.html#a-quick-introduction](https://spack.readthedocs.io/en/latest/containers.html#a-quick-introduction) |
| **ULHPC Docs** | ULHPC-Docs Software Config | Reference to a pull request on the ULHPC documentation for software stack configuration. | [https://github.com/ULHPC/ulhpc-docs/pull/136/files#diff-b8c99d3043663cc8a810ccef8174e284919267300179122a03e0526ad2d30a8aR128](https://github.com/ULHPC/ulhpc-docs/pull/136/files#diff-b8c99d3043663cc8a810ccef8174e284919267300179122a03e0526ad2d30a8aR128) |
| **ULHPC Spack Config** | Spack vs. ULHPC Std Env | GitHub Gist comparing Spack configuration against the standard software set on uni.lu clusters. | [https://gist.github.com/jhale/23d4d7646e2dc05d0adc0395767d044a](https://gist.github.com/jhale/23d4d7646e2dc05d0adc0395767d044a) |
| **GIT** | GIT Documentation | Official documentation for version control. | [https://docs.github.com/en/get-started](https://docs.github.com/en/get-started) |
| **GIT** | GIT Reference | Comprehensive command line reference. | [https://git-scm.com/docs](https://git-scm.com/docs) |
| **GIT** | GIT Interactive Tutorial | Learn Git branching through an interactive web tool. | [https://learngitbranching.js.org/?locale=fr\_FR](https://learngitbranching.js.org/?locale=fr_FR) |
| **FE Library** | Deal.ii Website | Official website for the Deal.ii finite element library. | [https://www.dealii.org/](https://www.dealii.org/) |
| **FE Library** | Deal.ii Tutorials | Collection of official tutorials. | [https://www.dealii.org/current/doxygen/deal.II/Tutorial.html](https://www.dealii.org/current/doxygen/deal.II/Tutorial.html) |
| **FE Library** | Deal.ii Documentation | Main documentation page. | [https://www.dealii.org/current/index.html](https://www.dealii.org/current/index.html) |
| **FE Library** | Deal.ii Library Reference | Doxygen reference for the C++ library. | [https://www.dealii.org/current/doxygen/deal.II/index.html](https://www.dealii.org/current/doxygen/deal.II/index.html) |
| **Meshing Tool** | Pygmsh | Python library for easy GMSH scripting. | [https://pypi.org/project/pygmsh/](https://pypi.org/project/pygmsh/) |
| **Meshing Tool** | Pygmsh GitHub | Source code repository. | [https://github.com/nschloe/pygmsh](https://github.com/nschloe/pygmsh) |
| **Meshing Tool** | Meshio | Tool for reading and writing various mesh formats. | [https://pypi.org/project/meshio/](https://pypi.org/project/meshio/) |
| **Meshing Tool** | Meshio GitHub | Source code repository. | [https://github.com/nschloe/meshio](https://github.com/nschloe/meshio) |
| **Meshing Tool** | Pymesh | Geometry processing library for Python. | [https://pymesh.readthedocs.io/en/latest/](https://pymesh.readthedocs.io/en/latest/) |
| **Meshing Tool** | Pymesh GitHub | Source code repository. | [https://github.com/PyMesh/PyMesh](https://github.com/PyMesh/PyMesh) |
| **Meshing Tool** | Netgen/NGSolve | Finite element code and mesher. | [https://ngsolve.org/](https://ngsolve.org/) |
| **Meshing Tool** | NGSolve Tutorials | Tutorials by J. Schöberl. | [https://jschoeberl.github.io/iFEM/intro.html](https://jschoeberl.github.io/iFEM/intro.html) |
| **Meshing Tool** | Neper | Tool for generating polycrystals/complex microstructures. | [https://neper.info/](https://neper.info/) |
| **Boolean Lib** | CGAL | Computational Geometry Algorithms Library. | [https://github.com/CGAL/cgal](https://github.com/CGAL/cgal) |
| **Boolean Lib** | VTK-based Booleans | Improved boolean operations using VTK. | [https://github.com/zippy84/vtkbool](https://github.com/zippy84/vtkbool) |
| **Boolean Lib** | MeshLib | Mesh processing library. | [https://github.com/MeshInspector/MeshLib](https://github.com/MeshInspector/MeshLib) |
| **Boolean Lib** | mcut | Mesh cutting utility. | [https://github.com/cutdigital/mcut](https://github.com/cutdigital/mcut) |
| **Boolean Lib** | trimesh | Python library for loading and using meshes. | [https://github.com/mikedh/trimesh](https://github.com/mikedh/trimesh) |
| **Boolean Lib** | vedo (Boolean Example) | Python library for 3D visualization and analysis. | [https://github.com/marcomusy/vedo/blob/master/examples/basic/boolean.py](https://github.com/marcomusy/vedo/blob/master/examples/basic/boolean.py) |
| **Boolean Lib** | PyMeshLab | Python interface to the MeshLab meshing system. | [https://github.com/cnr-isti-vclab/PyMeshLab](https://github.com/cnr-isti-vclab/PyMeshLab) |
| **Boolean Lib** | Pycork (PyPI) | Python bindings for the Cork boolean library. | [https://pypi.org/project/pycork/](https://pypi.org/project/pycork/) |
| **Boolean Lib** | Cork (C Library) | The underlying C library for Pycork. | [https://github.com/gilbo/cork](https://github.com/gilbo/cork) |
| **Boolean Lib** | Blender | Open-source 3D creation suite (includes Python module). | [https://www.blender.org/](https://www.blender.org/) |
| **Tomography** | Tomo2FE GitHub | Workflow for converting tomography data to compliant multipart meshes. | [https://github.com/ANR-MultiFIRE/TomoToFE/blob/main/workflow2/Workflow2-Python.ipynb](https://github.com/ANR-MultiFIRE/TomoToFE/blob/main/workflow2/Workflow2-Python.ipynb) |
| **Tomography** | Tomo2FE Article | Associated RILEM technical article. | [https://letters.rilem.net/index.php/rilem/article/view/184](https://letters.rilem.net/index.php/rilem/article/view/184) |
| **Clean Code** | Clean Code - Lesson 1 | Robert C. Martin's "Uncle Bob" series on clean code principles. | [https://www.youtube.com/watch?v=7EmboKQH8lM](https://www.youtube.com/watch?v=7EmboKQH8lM) |
| **Clean Code** | Clean Code - Conference Series | Playlist of related Clean Code conference talks. | [https://www.youtube.com/watch?v=7EmboKQH8lM\&list=PLmmYSbUCWJ4x1GO839azG\_BBw8rkh-zOj\&index=1](https://www.youtube.com/watch?v=7EmboKQH8lM&list=PLmmYSbUCWJ4x1GO839azG_BBw8rkh-zOj&index=1) |
| **Documentation**| Doxygen | Tool for generating documentation from source code. | [https://www.doxygen.nl/index.html](https://www.doxygen.nl/index.html) |



### General Coding Remark

If you are unable to find the documentation for a specific object, you can use the interactive `help()` command within **ipython3**. For instance, if you have a class object named `mesh`, typing `mesh.` and pressing **Tab** will allow you to explore its attributes, after which you can execute `help(mesh.attribute)` to get its specific documentation.

### Ubuntu Note
On **Ubuntu jammy**, a known conflict with FEniCSx v0.8.0 may require the removal of the `python3-numba` package using the following command:
```sh
sudo apt remove python3-numba

```
