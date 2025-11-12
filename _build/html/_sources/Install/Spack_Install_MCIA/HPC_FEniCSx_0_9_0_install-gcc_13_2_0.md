# 📖 FEniCSx v0.9.0 Installation on HPC (MCIA)

This document provides a comprehensive, step-by-step guide for installing **FEniCSx v0.9.0** using Spack on an HPC cluster, specifically referencing the **MCIA** environment.

This procedure utilizes external system modules for the compiler (`gcc@13.2.0`) and build tools (`cmake@3.27.6`), which is common practice on shared clusters.

> **⚠️ Resource Warning:** Installing the complete environment (including `gmsh` and `opencascade`) may require significant resources (e.g., **$\ge$ 16 GB RAM**) and can last for **more than 8 hours**, particularly if dependencies like LLVM need to be recompiled. A minimal `spack.yaml` is also provided for quick test.

## ⚡ Quick Use Guide (After Successful Installation)
This section details how to load and use the FEniCSx environment, plus handy connection aliases for the cluster.

### 🌐 Connect to the HPC 
You can create an alias to place in your local machine's ~/.bash_aliases (or equivalent) to simplify connection to the cluster.
```bash
# Connect to the HPC cluster
alias hpc_connect='ssh <username>@curta.mcia.fr'
```

>*Note:* Remember to replace <username> with your actual cluster username.

### 📁 Access files online (upload/download)
The [https://curta3.mcia.fr/](https://curta3.mcia.fr/) web service allows users to access Curta from a web browser using the Open OnDemand software.


### 📁 Optional: Local File System Mount (Requires sshfs)
These aliases use sshfs (SSH Filesystem) to mount your remote directories directly to a local folder (e.g., `$HOME/curta`), allowing you to browse and manage files as if they were local.

>*Prerequisite:* You must have the sshfs tool installed on your local computer.
```bash
# Mount the remote HOME directory to $HOME/curta
alias hpc_home_mount='sshfs <username>@curta.mcia.fr:/gpfs/home/<username> $HOME/curta'

# Mount the remote SCRATCH directory to $HOME/curta
alias hpc_scratch_mount='sshfs <username>@curta.mcia.fr:/scratch/<username> $HOME/curta'

# Unmount the local folder
alias hpc_unmount='fusermount -u $HOME/curta'
```

>*Tip:* Before running a mount alias, ensure the local directory (`$HOME/curta`) exists: `mkdir -p $HOME/curta`.

### 💻 Running a Job: Sbatch File Template

This template demonstrates how to properly configure a Slurm batch script to load the installed FEniCSx module and execute a parallel simulation.
```bash
#!/bin/bash

# ==============================================================================
# 1. SLURM DIRECTIVES (MODIFY THESE FOR YOUR JOB)
# ==============================================================================
#SBATCH --job-name=fenicsx_test_two_nodes
#SBATCH --constraint=compute
#SBATCH --nodes=2
#SBATCH --ntasks=4                      # Total number of MPI processes (4 total)
#SBATCH --tasks-per-node=2              # Processes per node (2 per node on 2 nodes = 4 total)
#SBATCH --time=0-00:35:00               # Short runtime
#SBATCH --mem=8GB                       # Total memory reserved for the job
##SBATCH --mem-per-cpu=8g              # Alternative: Request memory per core
#SBATCH --output=slurm_test_%j.out
#SBATCH --error=slurm_test_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail>@<domain>.fr # Replace with your email address
#SBATCH --chdir=/gpfs/home/<USERNAME>/your_working_dir # Set the working directory

# Exit immediately if a command exits with a non-zero status
set -e

# useful informations to print
echo "==============================================================================" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "==============================================================================" 

echo "--- STARTING FENICSX MODULE TEST ---"

# ==============================================================================
# 2. MODULE LOADING
# ==============================================================================
echo "Loading FEniCSx module (includes Python and Intel MPI)..."
# Ensure a clean environment before loading the FEniCSx module
module purge

# Load the installed FEniCSx environment module
# This command automatically sets all PATHs, LD_LIBRARY_PATHs, and Python paths.
module load fenicsx/0.9.0

# ==============================================================================
# 3. EXECUTION 
# ==============================================================================

# Print job parameters for verification
echo "Starting execution on $(hostname) at $(date)"
echo "Number of nodes requested: $SLURM_NNODES"
echo "Total MPI tasks: $SLURM_NTASKS"

echo "Elementary test..."
mpirun -n $SLURM_NTASKS python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}')"

# Execute the parallel Python script using mpirun
# mpirun uses the OpenMPI instance loaded by the FEniCSx module
echo "Running: mpirun -n $SLURM_NTASKS python <filename>.py"

# IMPORTANT: Ensure <filename>.py is in the directory specified by --chdir
mpirun -n $SLURM_NTASKS python <filename>.py

# ==============================================================================
# 4. POST-PROCESSING (OPTIONAL)
# ==============================================================================
echo "Job finished at $(date)"
```

### Using the Generated Module (Recommended)

1.  **Load the FEniCSx module** (after module propagation):
    ```bash
    module load fenicsx/0.9.0
    ```
2.  **Execute Parallel Job:** Use `mpirun` or `srun` within your Slurm script.

>*Tip:* Some test files are available at [Tutorials for test](https://github.com/Th0masLavigne/FEniCSx_GMSH_tutorials).

### Using Spack Environment Directly

This approach is required before the module file is deployed system-wide.

1.  **Source the Spack setup:**
    ```bash
    source /gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025/share/spack/setup-env.sh
    ```
2.  **Activate the environment:**
    ```bash
    spack env activate --dir /gpfs/softs/contrib/apps/fenicsx/spack-config-28.10.2025/envs/fenicsx/0.9.0/shared
    ```
3.  **Export for OpenMPI/PSM2 (Crucial on MCIA):**
    ```bash
    export OMPI_MCA_orte_precondition_transports=$(printf '%08x%08x-%08x%08x' \
        ${SLURM_JOB_ID:-0} ${SLURM_STEP_ID:-0} ${SLURM_JOB_ID:-0} ${SLURM_STEP_ID:-0})
    ```
4.  **Export for NLOPT Python Bindings (If installed):**
    ```bash
    export PYTHONPATH=/gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025/opt/spack/linux-skylake_avx512/nlopt-2.9.1-joq7eb4qexm46mtfyhukmc6zxabfcrve/lib64/python3.11/site-packages:$PYTHONPATH
    ```

---

## 1. Prerequisites and Interactive Session Setup

It is **highly recommended** to perform this installation within a long-running, stable interactive session.

### Starting an Interactive Session

Request sufficient resources for the build process (32GB RAM is recommended).

```bash
srun --nodes=1 --ntasks=16 --time=2-23:00:00 --mem=32G --pty /bin/bash -i
```

> **Tip:** If you believe your job has frozen during a long-running step, you can reconnect to the same node using job overlap: `srun --pty --overlap --jobid <jobid> /bin/bash -i`
> Find your `<jobid>` with `squeue -u $USER`. You can then monitor the node with `top -u $USER`.

Please note that the module deployment file is the only one that explcitly requires **user modifications** for paths and variables in case of changed version.

### Cleaning the Environment

To prevent conflicts with existing modules, especially OpenMPI versions, it is recommended to purge all loaded modules:

```bash
module purge
```

-----

## 2\. Step 1: Define Environment Variables and Logging 📝

We define all key installation paths and set up logging to capture almost all terminal output.

> **Note:** The use of `2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}` is appended to most commands below to ensure most output (stdout and stderr) are well captured in the log file, while also displaying it in the console. However, `export` and `source` does not work well if it is added.

```bash
# ==================================================
# Create the export variables
# ==================================================

# Specify the name of your log file.
export REPORT_LOG_FILENAME="report_install_fenicsx_$(date +'%Y-%m-%d_%H-%M-%S')_job-id_${SLURM_JOBID}.md"
export REPORT_LOG_DIRECTORY="/gpfs/home/tlavigne002"
echo "Starting log file creation at: ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}" | tee ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Choose the FEniCSx version to install
export FENICSX_VERSION="0.9.0"

# Path for cloning the Spack repository (including packages.yaml)
export FENICSX_SPACK_ROOT="/gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025"

# Path to packages.yaml defining externals
export SPACK_PACKAGES_YAML="${FENICSX_SPACK_ROOT}/etc/spack/packages.yaml"

# Root path for the Spack environment configuration
export FENICSX_SPACK_CONFIG_ROOT="/gpfs/softs/contrib/apps/fenicsx/spack-config-28.10.2025"

# Directory where the Spack environment will be created (contains spack.yaml)
export SPACK_ENV_DIR="${FENICSX_SPACK_CONFIG_ROOT}/envs/fenicsx/${FENICSX_VERSION}/shared"

# Prefix for where the Tcl/Lmod modulefiles will be deployed
export MODULEFILE_PREFIX="/gpfs/softs/contrib/modulefiles/fenicsx"

# Use the number of cores allocated by slurm for parallel compilation
export NCORES=$SLURM_NTASKS

echo "All variables set." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

-----

## 3\. Step 2: Configure External Modules (GCC and CMake)

We need to identify the exact install paths for our external modules to correctly define them in `packages.yaml`.

On the cluster you can check the available modules using `module avail`. Based on this exhaustive list you can identify the external packages you want to work with but be careful of potential conflicts if different compilers were used, or about consistency between verions.

In our case, we consider only gcc@13.2.0 and cmake@3.27.6. You can either automate the path identification to limit user error or copy paste them directly in tour packages.yaml by checking `module show <name_of_the_module>`.

> *Note:* For Openmpi, we use psm2 communication here, with pmix.

### Automated Path Detection 🔗
We use a helper function to reliably identify and export the paths for the external `gcc@13.2.0` and `cmake@3.27.6` modules. Execute the following block to find the exact prefix paths for the external modules:

```bash
get_module_path() {
    module show "$1" 2>&1 \
        | grep -E '(/gpfs/softs/[^[:space:]]+)' -o \
        | grep -v '/modulefiles/' \
        | head -1 \
        | sed 's/\/bin$//' \
        | sed 's/:$//'
}

# Identify paths for CMake and GCC
export PREFIX_CMAKE=$(get_module_path cmake/3.27.6)
export PREFIX_GCC=$(get_module_path compilers/gcc/13.2.0)

echo "CMAKE prefix identified: ${PREFIX_CMAKE}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
echo "GCC prefix identified: ${PREFIX_GCC}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Ensure clean workspace
module purge
```

### Forcing CMake Path 
(Crucial Fix, see: [Spack Installation Issue for FEniCSx 0.9.0 (py-fenics-basix): CMake Version Mismatch - installation - FEniCS Project](https://fenicsproject.discourse.group/t/spack-installation-issue-for-fenicsx-0-9-0-py-fenics-basix-cmake-version-mismatch/18330))

This step resolves a common Spack/Python build issue where a system-wide CMake is chosen instead of the correct module version.

```bash
# Ensure force SKBUILD_CMAKE_EXECUTABLE for Python packages
export SKBUILD_CMAKE_EXECUTABLE="${PREFIX_CMAKE}/bin/cmake"
echo export SKBUILD_CMAKE_EXECUTABLE="${SKBUILD_CMAKE_EXECUTABLE}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Ensure force CMAKE_EXECUTABLE for Spack itself
export CMAKE_EXECUTABLE="${PREFIX_CMAKE}/bin/cmake"
echo export CMAKE_EXECUTABLE="${CMAKE_EXECUTABLE}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Ensure the correct CMake is at the front of the PATH
export PATH="${PREFIX_CMAKE}/bin:${PATH}"
echo export PATH="${PATH}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

-----

## 4\. Step 3: Initialize and Configure Spack ⚙️

### Spack Setup and Initialization

Clone Spack and source the environment setup script.

```bash
# Create necessary directories
mkdir -p ${FENICSX_SPACK_ROOT} 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
mkdir -p ${FENICSX_SPACK_CONFIG_ROOT} 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
mkdir -p ${SPACK_ENV_DIR} 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Clone or Update Spack
if [ ! -d "${FENICSX_SPACK_ROOT}/.git" ]; then
    echo "Cloning Spack (latest from git) into ${FENICSX_SPACK_ROOT}..." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
    git clone --depth=2 https://github.com/spack/spack.git ${FENICSX_SPACK_ROOT} 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
else
    echo "Spack is present. Updating..." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
    (cd ${FENICSX_SPACK_ROOT} && git pull) 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
fi

# Disable local config to limit conflicts
export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CONFIG_PATH=${FENICSX_SPACK_CONFIG_ROOT}

# Enable spack shell
source ${FENICSX_SPACK_ROOT}/share/spack/setup-env.sh 

echo "Spack initialized (version: $(spack --version))" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
spack debug report 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

### Create `packages.yaml` for Externals

We define the externally provided compiler (`gcc`), build system (`cmake`), and job scheduler (`slurm`). Note the explicit paths for the compiler executables (`c`, `cxx`, `fortran`) and the mandatory `languages='c,c++,fortran'` for the compiler.

```bash
echo "Creating packages.yaml at ${SPACK_PACKAGES_YAML}..." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

cat << EOF > ${SPACK_PACKAGES_YAML} 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
packages:
  # 1. Compiler (keep modules as it is a compiler)
  gcc:
    externals:
    # include binutils. Beware, according to the ressources, need for binutils >=2.40
    - spec: gcc@13.2.0+binutils languages='c,c++,fortran'
      # From our variable
      prefix: ${PREFIX_GCC}
      extra_attributes:
        compilers:
          c: ${PREFIX_GCC}/bin/gcc
          cxx: ${PREFIX_GCC}/bin/g++
          fortran: ${PREFIX_GCC}/bin/gfortran
      modules:
      - compilers/gcc/13.2.0
    buildable: false
  cmake:
    externals:
    - spec: cmake@3.27.6
      prefix: ${PREFIX_CMAKE}
      modules: 
      - cmake/3.27.6
    buildable: false
  slurm:
    externals:
    - spec: slurm@23.11.5 sysconfdir=/etc/slurm
      prefix: /usr
    buildable: false
EOF

echo "packages.yaml created at ${SPACK_PACKAGES_YAML}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
echo "Content of packages.yaml:" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
cat ${SPACK_PACKAGES_YAML} 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

-----

## 5\. Step 4: Define the FEniCSx Spack Environment (`spack.yaml`) 🌿

This defines the full dependency graph for the complete FEniCSx stack, including visualization and extra physics libraries required by our team.

### Complete `spack.yaml` Configuration

This is the complete environment set for MCIA, including specific librairies required by some PhDs of the team:

```bash
echo "Creating complete spack.yaml at ${SPACK_ENV_DIR}/spack.yaml..." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

cat << EOF > ${SPACK_ENV_DIR}/spack.yaml 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
# spack.yaml for fenicsx/${FENICSX_VERSION}/shared
spack:
  specs:
  # FEniCSx Core Components
  - fenics-dolfinx@${FENICSX_VERSION}+adios2+petsc
  - py-fenics-dolfinx@${FENICSX_VERSION}+petsc4py+slepc4py

  - py-fenics-ffcx@${FENICSX_VERSION}
  - fenics-ufcx
  - py-fenics-ufl@2024.2.0
  
  - py-fenics-basix@${FENICSX_VERSION}
  - fenics-basix@${FENICSX_VERSION}

  # To respect omnipath, add fabrics psm2 and slurm management and PMI
  - openmpi@5.0.8 fabrics=psm2 schedulers=slurm 
  
  # Other Dependencies
  - petsc+mumps+fortran+superlu-dist~trilinos
  - adios2~sst+python
  - gmsh+opencascade
  - vtk@9.5.1+mpi+python+shared
  - hdf5+hl+shared+mpi #by default, hdf5 was without hl causing crash
  
  # Common Python Libraries
  - python
  - py-pip
  - py-setuptools
  - py-wheel
  - python-venv
  - py-gmsh
  - py-pygmsh
  - py-numba
  - py-scipy
  - py-statsmodels
  - py-matplotlib
  - py-imageio
  - py-adios4dolfinx
  - py-trimesh
  - py-rtree
  - py-numpy
  - py-pandas
  - py-pyvista
  - py-h5py

  # Other Libraries
  - nlopt+python
  - neper
  - py-pytz
  - py-scikit-optimize
  - py-seaborn
  - py-python-pptx

  config:
    concretizer:
      unify: true
    env_vars:
      # Crucial: Override environment variables to force the correct CMake path
      set:
        SKBUILD_CMAKE_EXECUTABLE: ${PREFIX_CMAKE}/bin/cmake
        CMAKE_EXECUTABLE: ${PREFIX_CMAKE}/bin/cmake
        PATH: ${PREFIX_CMAKE}/bin:$PATH
EOF

echo "Content of spack.yaml:" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
cat ${SPACK_ENV_DIR}/spack.yaml 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

> *Note:* You can browse [Spack Packages](https://packages.spack.io/) to discover additional software available for installation.

> **Maintenance Note:** If you modify `packages.yaml` or `spack.yaml`, you must run:
>
> ```bash
> spack env deactivate
> spack env activate -p --dir ${SPACK_ENV_DIR}
> spack deconcretize -a -y
> spack concretize
> ```
> By doing so, you ensure having the updated spack tree.

-----

## 6\. Step 5: Install the environment

We activate the environment, resolve dependencies (concretize), and start the installation.

```bash
# Activate the environment. Use -p for persistent activation (recommended for interactive shells).
echo "spack env activate -p --dir ${SPACK_ENV_DIR}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
spack env activate -p --dir ${SPACK_ENV_DIR}

# Concretize the environment to resolve the full dependency graph
echo "Concretizing environment..." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
spack concretize 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Verify the externals (GCC/CMake) were respected
echo "Verify the externals are respected (should show 'E' for externals):" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
spack find -c 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Verify dependencies for a critical component (e.g., py-fenics-basix)
echo "Verify the cmake for basix:" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
spack spec py-fenics-basix 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Start the installation
echo "Starting FEniCSx installation (using ${NCORES} cores for -j)..." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
# Replace %gcc@13.2.0 by yours if needed.
spack install -j${NCORES} %gcc@13.2.0 --fail-fast 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

echo "FEniCSx installation complete. $(date +'%Y-%m-%d_%H-%M-%S'). View installed packages at: ${FENICSX_SPACK_CONFIG_ROOT}/envs/fenicsx/0.9.0/shared/.spack-env/view" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

> *Note:* For VTK or LLVM for example, the installation times are quite long (up to 12 hours each). 


If you need more info about where is the package installed to check for example if multiply installed, run:
```bash
spack location -i openmpi
spack find -p openmpi
```

> **Rebuilding:** If you need to rebuild core FEniCSx components from source (e.g., after updating Spack), use the `--overwrite` flag:
> `spack install --overwrite -j${NCORES} fenics-basix py-fenics-basix py-fenics-ffcx fenics-ufcx py-fenics-ufl fenics-dolfinx py-fenics-dolfinx`


Some packages might fail to install probably due to fetching step. For these missing packages identified with `spack find -c`, you can reinstall immediatly them using (and now everything should be installed). For example: 
```bash
spack install -j${NCORES} %gcc@13.2.0 gsl@2.8 neper@4.10.1 py-adios4dolfinx@0.9.4 py-certifi@2025.7.14 py-charset-normalizer@3.4.4 py-gmsh@4.13.1 py-h5py@3.13.0 py-idna@3.10 py-imageio@2.37.0 py-mako@1.3.10 py-meshio@5.0.1 py-pkgconfig@1.5.5 py-platformdirs@4.4.0 py-poetry-core@2.2.0 py-pooch@1.8.2 py-pyaml@21.8.3 py-pygmsh@7.1.17 py-python-pptx@0.6.23 py-pyvista@0.46.3 py-pyyaml@6.0.3 py-requests@2.32.5 py-rtree@1.4.1 py-scikit-optimize@0.9.0 py-scooby@0.10.0 py-seaborn@0.13.2 py-trimesh@3.17.1 py-typing-extensions@4.14.1 py-urllib3@2.5.0 py-xlsxwriter@3.1.7 xcb-proto@1.17.0
```

> **Handling Installation Failures:** If packages fail (often due to timeouts/network issues), re-force installation for the specific failed packages.
>
> ```bash
> # Example to re-force installation for failed Python packages:
> spack install -j${NCORES} %gcc@13.2.0 py-requests@2.32.5 py-idna@3.10 py-rtree@1.4.1 [...]
> ```



-----

## 7\. Step 6: Testing and Validation ✅

Once installed, confirm the environment is functional, paying close attention to MPI/PSM2 configuration.


### Basic MPI Test

Ensure you are using the right psm2 communication:
```bash
ompi_info --param mtl psm2 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
# expect MCA mtl: psm2 (MCA v2.1.0, API v2.0.0, Component v5.0.8)
```

This test verifies `dolfinx` and `mpi4py` are correctly linked and running in parallel. 

>Based on the following note: if you use spack directly and not module load, you need to define the following variable  before running mpirun (16 characters is required, see [workaround pmix with psm2](https://github.com/open-mpi/ompi/issues/13397).

```bash
# Set the crucial PSM2/PMIX key if the module is not yet loaded
if [[ -n "$SLURM_JOB_ID" && -n "$SLURM_STEP_ID" ]]; then
    export OMPI_MCA_orte_precondition_transports=$(printf '%08x%08x-%08x%08x' \
        ${SLURM_JOB_ID:-0} ${SLURM_STEP_ID:-0} ${SLURM_JOB_ID:-0} ${SLURM_STEP_ID:-0})
fi

mpirun -n $SLURM_NTASKS python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}')" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

### Debugging PMIX/PSM2 Errors

If the PMIX transport error persists, test with a fallback MCA configuration:

```bash
mpirun -n 4 --mca pml ob1 --mca btl self,vader --mca orte_precondition_transports false python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}')"  2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

### Missing Python bindings:

If a python binding is missing like for nlopt, if you use spack directly (and not module load that handles it already), add:
```bash
export PYTHONPATH=/gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025/opt/spack/linux-skylake_avx512/nlopt-2.9.1-joq7eb4qexm46mtfyhukmc6zxabfcrve/lib64/python3.11/site-packages:$PYTHONPATH
```

-----

## 8\. Step 7: Module Deployment (Propagating to Tcl on MCIA)

This step generates the final Tcl modulefile for system-wide access. This file must be placed in a directory sourced by the cluster's module system. Notably it sets the missing variables automatically and loads the view (i.e. the whole spack environment) under a single flag: `fenicsx/0.9.0`.


Create the target directory for the modulefile:
```bash
mkdir -p "${MODULEFILE_PREFIX}"
```

Create the TCL module file. The use of 'EOF' ensures that variables like ${FENICSX_VERSION} are expanded now, while internal logic like the 'if {[info exists env(SLURM_JOB_ID)]}' remains as literal text. This file is pointing to your spack install view and auto-generate the unique key for psm2/pmix. It also specifies some python bindings if missing like the one of nlopt in our case:
```bash
cat << 'EOF' > ${MODULEFILE_PREFIX}/${FENICSX_VERSION}
#%Module1.0#####################################################################
##
## FEniCSx 0.9.0 — modulefile (HPC shared Spack environment)
##

proc ModulesHelp { } {
    puts stderr "Loads FEniCSx ${FENICSX_VERSION} environment built with Spack"
    puts stderr "Includes Dolfinx, Basix, UFL, ADIOS2, PETSc, etc."
}
module-whatis "FEniCSx 0.9.0 — Shared HPC Spack environment (Dolfinx stack)"

# ============================================================
# Core environment variables
# ============================================================
setenv FENICSX_VERSION "0.9.0"
setenv FENICSX_ENV_ROOT "/gpfs/softs/contrib/apps/fenicsx/spack-config-28.10.2025/envs/fenicsx/0.9.0/shared"
setenv FENICSX_VIEW     "$env(FENICSX_ENV_ROOT)/.spack-env/view"

# ============================================================
# Path configuration
# ============================================================

prepend-path PATH              "$env(FENICSX_VIEW)/bin"
prepend-path LD_LIBRARY_PATH   "$env(FENICSX_VIEW)/lib"
prepend-path LD_LIBRARY_PATH   "$env(FENICSX_VIEW)/lib64"
prepend-path CMAKE_PREFIX_PATH "$env(FENICSX_VIEW)"
prepend-path PKG_CONFIG_PATH   "$env(FENICSX_VIEW)/lib/pkgconfig"
prepend-path PKG_CONFIG_PATH   "$env(FENICSX_VIEW)/lib64/pkgconfig"
prepend-path CPATH             "$env(FENICSX_VIEW)/include"
prepend-path MANPATH           "$env(FENICSX_VIEW)/share/man"

# Python integration
prepend-path PYTHONPATH "$env(FENICSX_VIEW)/lib/python3.11/site-packages"
#
# NLOpt Python bindings (manual install outside view)
prepend-path PYTHONPATH "/gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025/opt/spack/linux-skylake_avx512/nlopt-2.9.1-joq7eb4qexm46mtfyhukmc6zxabfcrve/lib64/python3.11/site-packages"

# ============================================================
# UCX / PSM2 runtime tuning
# ============================================================

# Recommended UCX transport settings for Omni-Path
# setenv UCX_TLS "psm2,sm,self"
# setenv UCX_NET_DEVICES "all"

# ============================================================
# Optional: ensure unique PSM2 key under Slurm because Error obtaining unique transport key from PMIX (OMPI_MCA_orte_precondition_transports not present in the environment). when loading fenicsx========================

if {[info exists env(SLURM_JOB_ID)]} {
    # Default to step ID 0, which is correct for the main sbatch script shell
    set stepid 0
    if {[info exists env(SLURM_STEP_ID)]} {
        set stepid $env(SLURM_STEP_ID)
    }
    
    set jobid $env(SLURM_JOB_ID)
    
    # Key generation using Job ID and Step ID (defaulted to 0)
    set key [format "%08x%08x-%08x%08x" $jobid $stepid $jobid $stepid]
    setenv OMPI_MCA_orte_precondition_transports $key
}

# ============================================================
# Environment info
# ============================================================

setenv FENICSX_HOME "$env(FENICSX_VIEW)"
setenv FENICSX_SPACK_ENV "$env(FENICSX_ENV_ROOT)"

puts stderr ">> Loaded FEniCSx 0.9.0 (shared Spack environment)"
puts stderr ">> Using view: $env(FENICSX_VIEW)"
puts stderr ">> Using OMPI_MCA_orte_precondition_transports: $key"

EOF
```

> Trick for sbatch stability: STEP_ID not always defined => need to default it if not existing otherwise no key generated.

-----

## 9\. Cleanup and Reinstallation (For Maintenance) 🗑️

These commands are for maintenance purposes if a complete fresh start is required. **Run these manually outside of the installation log script.**

The spack got broken because of a single package (let's say pip install XXX). Reinstall only python and bindings:
> \*Uninstalling a single package but its dependencies too for reinstall clean:  `spack uninstall  --dependents -a python`

Full remove:
```bash
# 1. Uninstall all packages in the environment
spack uninstall -j${NCORES} %gcc@13.2.0 --fail-fast 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# 2. Clean up temporary stages
spack clean -s 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
spack clean -a 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# 3. Deactivate the environment
spack env deactivate 

# 4. Remove directories
echo "Removing Spack environment: ${SPACK_ENV_DIR}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
rm -rf ${SPACK_ENV_DIR}

echo "Removing Spack configuration root: ${FENICSX_SPACK_CONFIG_ROOT}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
rm -rf ${FENICSX_SPACK_CONFIG_ROOT}

echo "Removing Spack source repository: ${FENICSX_SPACK_ROOT}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
rm -rf ${FENICSX_SPACK_ROOT}
```

Caution to remove your package only:
```bash
echo "Removing Modules: ${MODULEFILE_PREFIX}/${FENICSX_VERSION}" 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
rm -rf ${MODULEFILE_PREFIX}/${FENICSX_VERSION}
```

-----

## Resources

  * Spack Tutorial: [spack-tutorial.readthedocs.io](https://spack-tutorial.readthedocs.io/en/latest/tutorial_configuration.html)

  * Jack Hale Spack/Module Example: [gist.github.com/jhale/23d4d7646e2dc05d0adc0395767d044a](https://gist.github.com/jhale/23d4d7646e2dc05d0adc0395767d044a)

  * FEniCSx Tutorials for Testing: [github.com/Th0masLavigne/FEniCSx\_GMSH\_tutorials](https://github.com/Th0masLavigne/FEniCSx_GMSH_tutorials)

  * FEniCSx Discourse (CMake issue discussion): [fenicsproject.discourse.group](https://www.google.com/search?q=https://fenicsproject.discourse.group/t/spack-installation-issue-for-fenicsx-0-9-0-py-fenics-basix-cmake-version-mismatch/18330/2)

  * Spack Packages: [packages.spack.io](https://packages.spack.io/)

  * Set variable for openmpi: [openmpi/pmix issue](https://github.com/open-mpi/ompi/issues/13397)



-----

#### Appendix: (Minimal) Configuration file
If you wish to only test the core FEniCSx installation without `gmsh` or other optional libraries, you can use a minimal configuration that excludes `gmsh+opencascade`, `py-gmsh`, `py-numba`, etc.
```bash
cat << EOF > ${SPACK_ENV_DIR}/spack_minimal.yaml 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
# spack_minimal.yaml for fenicsx/${FENICSX_VERSION}/shared
spack:
  specs:
  - python
  - py-pip
  - py-setuptools
  - py-wheel
  - python-venv
  - fenics-dolfinx@${FENICSX_VERSION}+adios2+petsc
  - py-fenics-dolfinx@${FENICSX_VERSION}+petsc4py+slepc4py
  - py-fenics-basix@${FENICSX_VERSION}
  - fenics-basix@${FENICSX_VERSION}
  - py-fenics-ufl@2024.2.0
  - py-fenics-ffcx@${FENICSX_VERSION}
  - petsc+mumps+fortran+superlu-dist~trilinos
  - adios2~sst+python
  # To respect omnipath, add fabrics psm2
  - openmpi@5.0.8 fabrics=psm2 schedulers=slurm

  config:
    concretizer:
      unify: true
    env_vars:
      set:
        SKBUILD_CMAKE_EXECUTABLE: ${PREFIX_CMAKE}/bin/cmake
        CMAKE_EXECUTABLE: ${PREFIX_CMAKE}/bin/cmake
        PATH: ${PREFIX_CMAKE}/bin:$PATH
EOF
```