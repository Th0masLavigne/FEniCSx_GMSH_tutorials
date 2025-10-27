# Installation Guide for FEniCSx v0.9.0

This guide provides detailed, copy-paste ready instructions to install **FEniCSx v0.9.0** and its dependencies across various operating systems, including dedicated automated workflows for **Windows (via WSL 2)** and **HPC environments (via Spack)**.

---

## 1. Linux Installation (Native or WSL2 Terminal)

For native Linux systems or when working directly inside a WSL2 distribution's terminal, you have two primary options: the APT repository for standard users or Spack for advanced/HPC users.

### 1.1. Local Installation of FEniCSx v0.9.0 (APT Method)

> **Note:** These instructions are tested for **Ubuntu 24.04 LTS**. This method installs FEniCSx directly into your system's packages.

#### Step 1 — Add the FEniCSx Repository
Access the package and update your shell environment:

```sh
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update
````

#### Step 2 — Install FEniCSx and Dependencies

Install FEniCSx at the specific version (`0.9.0.1`) and all necessary packages for meshing, plotting, and environment management.

```sh
# Install FEniCSx v0.9.0
sudo apt install -y 'fenicsx=2:0.9.0.1*'

# Install Python tools, Pandas, Gmsh, and rendering libraries for GUI/plotting
sudo apt install python3-pip python3-venv python3-pandas -y
sudo apt install gmsh xvfb libgl1-mesa-dri libglx-mesa0 libglu1-mesa mesa-utils -y
```

-----

## 2\. Windows Installation (Automated Batch Scripts)

For Windows users, using the provided batch scripts is the most reliable way to set up an isolated environment with FEniCSx. **Ensure Docker Desktop is installed and running for all Docker methods.**

### 2.1. Option A: WSL2 Docker Installer (Recommended for Consistency)

This approach uses **Docker within a dedicated WSL2 distribution** (`FEniCSxenv`), offering high performance and consistency.

1.  **Files:** Use the files from the **`WSL2_Docker_installer`** folder (containing `install_wsl2_fenicsx_docker.bat` and others).
2.  **Installation:** Run the installer from Windows (Command Prompt/PowerShell):
    ```bash
    install_wsl2_fenicsx_docker.bat
    ```
3.  **Usage:** After installation, launch your environment directly from Windows using the provided launchers:
    ```bash
    # Interactive shell (will prompt for a directory to mount)
    FEniCSx_interactive.bat

    # Run a command/script (edit the .bat file first to set script path)
    FEniCSx_command.bat
    ```
    *Alternatively, open your WSL terminal (`wsl -d FEniCSxenv`) and use the custom aliases:*
    ```bash
    fenicsx                       # Start interactive container shell
    fenicsx-run "python3 script.py" # Run non-interactive command
    ```

### 2.2. Option B: WSL2 Local Installer (Native `venv` in WSL)

This sets up a dedicated WSL distribution (`FEniCSxenv`) with a FEniCSx Python Virtual Environment (`venv`) that auto-activates upon login.

1.  **Files:** Use the files from the **`WSL_Local_installer`** folder (containing `install_local_fenicsx.bat`, `apt-installed-list.txt`, etc.).
2.  **Installation:** Run the installer from Windows (Command Prompt/PowerShell):
    ```bash
    install_local_fenicsx.bat
    ```
    *You will be prompted to enter the WSL username you created during the initial Ubuntu setup.*
3.  **Usage:** Launch the dedicated environment from Windows:
    ```bash
    wsl -d FEniCSxenv
    ```
    The `fenicsx-env` will be automatically active. Run scripts with:
    ```bash
    python3 filename.py
    ```

### 2.3. Option C: Windows Direct Docker

This method runs Docker commands directly from Windows, relying on Docker Desktop's native file mapping.

1.  **Prerequisites:** Ensure **Docker Desktop** is running.
2.  **Files:** Use the files from the **`Windows_docker_installer`** folder (containing `FEniCSx_docker_0-9-0.bat`, etc.).
3.  **Image Build:** Execute the build script. This will also integrate any custom packages found in the `src` folder.
    ```bash
    FEniCSx_docker_0-9-0.bat
    ```
4.  **Usage:**
    ```bash
    # Interactive shell
    Interactive_FEniCSx_0-9-0.bat

    # Non-interactive script run (e.g., Test.py or Main.py)
    FEniCSx_AGGREGATE.bat
    ```

-----

## 3\. HPC / Advanced Installation (Spack Package Manager)

Spack is the recommended tool for managing complex dependencies on High-Performance Computing (HPC) clusters, allowing for precise control over compilers, MPI, and linear algebra libraries like PETSc.

### Spack Installation Procedure (Detailed)

This procedure was tested on **Red Hat Enterprise Linux 8.10**.

#### Step 1 — Get and Configure Spack

Download and set up the Spack environment.

```sh
# 1. Get Spack v0.23.0
wget [https://github.com/spack/spack/releases/download/v0.23.0/spack-0.23.0.tar.gz](https://github.com/spack/spack/releases/download/v0.23.0/spack-0.23.0.tar.gz)

# 2. Untar the archive
tar -xf spack-0.23.0.tar.gz

# 3. Optional: Disable local config (recommended for multi-user/HPC setups)
export SPACK_DISABLE_LOCAL_CONFIG=true

# 4. Activate Spack
cd spack-0.23.0
source ./share/spack/setup-env.sh
```

#### Step 2 — Configure Compilers and Environment

Add your system compiler to Spack, and optionally install a newer GCC version if needed.

```sh
# 1. Add system compiler (must be pre-installed on the OS)
spack compiler find 

# 2. Optional: Install a modern compiler (e.g., GCC 12)
spack install gcc@12
spack compiler find $(spack location -i gcc@12)
```

#### Step 3 — Create and Install FEniCSx Environment

Create a dedicated Spack environment, add the required packages, and initiate the installation.

```sh
# 1. Create a Spack environment
spack env create MyFenicsxEnv

# 2. Activate the environment
spack env activate MyFenicsxEnv

# 3. Add Core FEniCSx dependencies and solvers (specify exact versions for reproducibility)
# PETSc config: Disable Trilinos, enable Fortran, Mumps, and SuperLU-Dist
spack add petsc@3.22.1 ~trilinos +fortran +mumps +superlu-dist
spack add python@3.11.9
spack add py-fenics-dolfinx@0.9.0
spack add py-mpi4py@4.0.1
spack add py-numpy@2.1.2
spack add py-matplotlib
spack add py-pandas
spack add py-gmsh

# 4. Install all packages in the environment
spack install
```

**Success\!** You have installed FEniCSx v0.9.0 and its parallel dependencies.

#### Step 4 — Usage and Activation

To use the environment in any new terminal session or submission script:

```sh
# 1. Activate Spack environment (add this to your .bashrc or .bash_profile)
source /path/to/spack-0.23.0/share/spack/setup-env.sh

# 2. Activate your custom environment
spack env activate MyFenicsxEnv

# 3. Run a simulation using MPI
mpirun -np 4 python3 my_fenicsx_script.py
```

> **HPC Batch Scripts (Slurm):** You can reference the file **`Slurm_spack_fenicsx_0_9_0.sh`** for an example of how to correctly source and activate the Spack environment within a Slurm job submission script.
