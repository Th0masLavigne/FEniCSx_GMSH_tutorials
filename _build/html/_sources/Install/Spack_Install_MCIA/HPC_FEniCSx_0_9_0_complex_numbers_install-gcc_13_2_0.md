# 📖 \[Complex number\] FEniCSx v0.9.0 Installation on HPC (MCIA)

This document provides a comprehensive, step-by-step guide for installing **FEniCSx v0.9.0** with **complex** number supports, using Spack on an HPC cluster, specifically referencing the **MCIA** environment. Conversely to the set-up provided for the real number support, this installation consider minimal FEniCSx environment (no nlopt, no gmsh, etc.).

As explained in [J. Dokken Complex Poisson Problem](https://jsdokken.com/dolfinx-tutorial/chapter1/complex_mode.html), the difference relies in the linear algebra that now needs to support complex numbers.

> **⚠️ Note:** Load fenicsx/0.9.0 or fenicsx_complex/0.9.0 depending on your need but the two environment do not coexist.

This procedure utilizes external system modules for the compiler (`gcc@13.2.0`) and build tools (`cmake@3.27.6`), which is common practice on shared clusters.

## ⚡ Quick Use Guide (After Successful Installation)
Please refer to [MCIA Installation FEnicsx 0.9.0 - Real numbers support](https://th0maslavigne.github.io/FEniCSx_GMSH_tutorials/Install/Spack_Install_MCIA/HPC_FEniCSx_0_9_0_install-gcc_13_2_0.html).

### 💻 Sbatch File Template

This template demonstrates how to properly configure a Slurm batch script to load the installed FEniCSx module and execute a parallel simulation.

>To run multiple simulations with a variable value like for the mpirun or a parameter of your code, see at the really bottom of the page.

```bash
#!/bin/bash

# ==============================================================================
# 1. SLURM DIRECTIVES (MODIFY THESE FOR YOUR JOB)
# ==============================================================================
#SBATCH --job-name=fenicsx_test_two_nodes
#SBATCH --constraint=compute
##SBATCH -p i2m,i2m-resources
#SBATCH --partition=i2m,i2m-resources
#SBATCH --nodes=1
#SBATCH --ntasks=4                      # Total number of MPI processes (4 total)
##SBATCH --tasks-per-node=2              # Processes per node (2 per node on 2 nodes = 4 total)
#SBATCH --time=0-01:00:00               # Short runtime
#SBATCH --mem=16GB                       # Total memory reserved for the job
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
module load fenicsx_complex/0.9.0

# ==============================================================================
# 3. EXECUTION 
# ==============================================================================

# Print job parameters for verification
echo "Starting execution on $(hostname) at $(date)"
echo "Number of nodes requested: $SLURM_NNODES"
echo "Total MPI tasks: $SLURM_NTASKS"

echo "Elementary test..."
mpirun -n $SLURM_NTASKS python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}');from petsc4py import PETSc; from dolfinx.fem.petsc import assemble_vector; print(PETSc.ScalarType);"

# Execute the parallel Python script using mpirun
# mpirun uses the OpenMPI instance loaded by the FEniCSx module
echo "Running: mpirun -n $SLURM_NTASKS python <filename>.py"

# IMPORTANT: Ensure <filename>.py is in the directory specified by --chdir
mpirun -n $SLURM_NTASKS python <filename>.py

# ==============================================================================
# 4. POST-PROCESSING (OPTIONAL)
# ==============================================================================
echo "Job finished at $(date)"

echo "After a job finishes you can see the real use of CPUs and memory"
echo "1. module load slurm/wrappers"
echo "2. seff <job-id>"
echo "3. module purge"
```

---

## Installation procedure followed for the MCIA cluster

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

### Step 1: Define Environment Variables and Logging 📝

We define all key installation paths and set up logging to capture almost all terminal output.

> **Note:** The use of `2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}` is appended to most commands below to ensure most output (stdout and stderr) are well captured in the log file, while also displaying it in the console. However, `export` and `source` does not work well if it is added.

```bash
# ==================================================
# Create the export variables
# ==================================================

# Specify the name of your log file.
export REPORT_LOG_FILENAME="report_install_fenicsx_complex_$(date +'%Y-%m-%d_%H-%M-%S')_job-id_${SLURM_JOBID}.md"
export REPORT_LOG_DIRECTORY="/gpfs/home/tlavigne002"
echo "Starting log file creation at: ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}" | tee ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}

# Choose the FEniCSx version to install
export FENICSX_VERSION="0.9.0"

# Path for cloning the Spack repository (including packages.yaml)
export FENICSX_SPACK_ROOT="/gpfs/softs/contrib/apps/fenicsx_complex/spack-src-17.11.2025"

# Path to packages.yaml defining externals
export SPACK_PACKAGES_YAML="${FENICSX_SPACK_ROOT}/etc/spack/packages.yaml"

# Root path for the Spack environment configuration
export FENICSX_SPACK_CONFIG_ROOT="/gpfs/softs/contrib/apps/fenicsx_complex/spack-config-17.11.2025"

# Directory where the Spack environment will be created (contains spack.yaml)
export SPACK_ENV_DIR="${FENICSX_SPACK_CONFIG_ROOT}/envs/fenicsx_complex/${FENICSX_VERSION}/shared"

# Prefix for where the Tcl/Lmod modulefiles will be deployed
export MODULEFILE_PREFIX="/gpfs/softs/contrib/modulefiles/fenicsx_complex"

# Use the number of cores allocated by slurm for parallel compilation
export NCORES=$SLURM_NTASKS

echo "All variables set." 2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

-----

### Step 2: Configure External Modules (GCC and CMake)

We need to identify the exact install paths for our external modules to correctly define them in `packages.yaml`.

On the cluster you can check the available modules using `module avail`. Based on this exhaustive list you can identify the external packages you want to work with but be careful of potential conflicts if different compilers were used, or about consistency between verions.

In our case, we consider only gcc@13.2.0 and cmake@3.27.6. You can either automate the path identification to limit user error or copy paste them directly in tour packages.yaml by checking `module show <name_of_the_module>`.

> *Note:* For Openmpi, we use psm2 communication here, with pmix.

#### Automated Path Detection 🔗
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

#### Forcing CMake Path 
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

### Step 3: Initialize and Configure Spack ⚙️

#### Spack Setup and Initialization

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

#### Create `packages.yaml` for Externals

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

### Step 4: Define the FEniCSx Spack Environment (`spack.yaml`) 🌿

This defines the full dependency graph for the complete FEniCSx stack, including visualization and extra physics libraries required by our team.


#### Minimal `spack.yaml` Configuration

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
  - petsc+mumps+fortran+superlu-dist~trilinos+complex
  - adios2~sst+python
  - hdf5+hl+shared+mpi #by default, hdf5 was without hl causing crash
  
  # Common Python Libraries
  - python
  - py-pip
  - py-setuptools
  - py-wheel
  - python-venv
  - py-scipy
  - py-matplotlib
  - py-numpy

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

### Step 5: Install the environment

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

If you need more info about where is the package installed to check for example if multiply installed, run:
```bash
spack location -i openmpi
spack find -p openmpi
```

> **Rebuilding:** If you need to rebuild core FEniCSx components from source (e.g., after updating Spack), use the `--overwrite` flag:
> `spack install --overwrite -j${NCORES} fenics-basix py-fenics-basix py-fenics-ffcx fenics-ufcx py-fenics-ufl fenics-dolfinx py-fenics-dolfinx`


Some packages might fail to install probably due to fetching step. For these missing packages identified with `spack find -c`, you can reinstall immediatly them using (and now everything should be installed). For example: 
```bash
spack install -j${NCORES} %gcc@13.2.0 gsl@2.8 py-gmsh@4.13.1 py-h5py@3.13.0 
```

> **Handling Installation Failures:** If packages fail (often due to timeouts/network issues), re-force installation for the specific failed packages.
>
> ```bash
> # Example to re-force installation for failed Python packages:
> spack install -j${NCORES} %gcc@13.2.0 py-requests@2.32.5 py-idna@3.10 py-rtree@1.4.1 [...]
> ```



-----

### Step 6: Testing and Validation ✅

Once installed, confirm the environment is functional, paying close attention to MPI/PSM2 configuration.


#### Basic MPI Test

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

#### Debugging PMIX/PSM2 Errors

If the PMIX transport error persists, test with a fallback MCA configuration:

```bash
mpirun -n 4 --mca pml ob1 --mca btl self,vader --mca orte_precondition_transports false python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}')"  2>&1 | tee -a ${REPORT_LOG_DIRECTORY}/${REPORT_LOG_FILENAME}
```

#### Missing Python bindings:

If a python binding is missing like for nlopt, if you use spack directly (and not module load that handles it already), add:
```bash
export PYTHONPATH=/gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025/opt/spack/linux-skylake_avx512/nlopt-2.9.1-joq7eb4qexm46mtfyhukmc6zxabfcrve/lib64/python3.11/site-packages:$PYTHONPATH
```

-----

### Step 7: Module Deployment (Propagating to Tcl on MCIA)

This step generates the final Tcl modulefile for system-wide access. This file must be placed in a directory sourced by the cluster's module system. Notably it sets the missing variables automatically and loads the view (i.e. the whole spack environment) under a single flag: `fenicsx_complex/0.9.0`.


Create the target directory for the modulefile:
```bash
mkdir -p "${MODULEFILE_PREFIX}"
```

Create the TCL module file. The use of 'EOF' ensures that variables like ${FENICSX_VERSION} are expanded now, while internal logic like the 'if {[info exists env(SLURM_JOB_ID)]}' remains as literal text. This file is pointing to your spack install view and auto-generate the unique key for psm2/pmix. It also specifies some python bindings if missing like the one of nlopt in our case:
```bash
cat << 'EOF' > ${MODULEFILE_PREFIX}/${FENICSX_VERSION}
#%Module1.0#####################################################################
##
## FEniCSx 0.9.0 complex number support— modulefile (HPC shared Spack environment)
##

proc ModulesHelp { } {
    puts stderr "Loads FEniCSx ${FENICSX_VERSION} complex number support environment built with Spack"
    puts stderr "Includes Dolfinx, Basix, UFL, ADIOS2, PETSc, etc."
}
module-whatis "FEniCSx 0.9.0 complex number support — Shared HPC Spack environment (Dolfinx stack)"

# ============================================================
# Core environment variables
# ============================================================
setenv FENICSX_VERSION "0.9.0"
setenv FENICSX_ENV_ROOT "/gpfs/softs/contrib/apps/fenicsx_complex/spack-config-17.11.2025/envs/fenicsx_complex/0.9.0/shared"
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
# prepend-path PYTHONPATH "/gpfs/softs/contrib/apps/fenicsx/spack-src-28.10.2025/opt/spack/linux-skylake_avx512/nlopt-2.9.1-joq7eb4qexm46mtfyhukmc6zxabfcrve/lib64/python3.11/site-packages"

# ============================================================
# UCX / PSM2 runtime tuning
# ============================================================

# Recommended UCX transport settings for Omni-Path
# setenv UCX_TLS "psm2,sm,self"
# setenv UCX_NET_DEVICES "all"

# ============================================================
# Optional: ensure unique PSM2 key under Slurm because Error obtaining unique transport key from PMIX (OMPI_MCA_orte_precondition_transports not present in the environment). when loading fenicsx_complex
# ============================================================

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

puts "echo '============================================================'"
puts "echo 'Loaded FEniCSx 0.9.0 complex number support'"
puts "echo '============================================================'"
puts "echo 'Using view: $env(FENICSX_VIEW)'"
puts "echo 'Using OMPI_MCA_orte_precondition_transports: $key'"
puts "echo  '============================================================'"
puts "echo  '============================================================'"

EOF
```

> Trick for sbatch stability: STEP_ID not always defined => need to default it if not existing otherwise no key generated.

-----

### Cleanup and Reinstallation (For Maintenance) 🗑️

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


## Array in slurm

This SLURM script uses the **Job Array** feature to efficiently run **800 independent, identical simulations** in parallel. Each simulation is an **MPI** (Message Passing Interface) task utilizing **16 CPU cores**. The job array automatically manages the submission and tracking of these 800 tasks (indexed **0 to 799**), making it ideal for high-throughput parametric studies (like FEniCSx simulations with varying parameters).

-----

### ⚙️ SLURM Array Configuration Explained

The core of the array configuration is defined by the following directives in section 1:

  * **`#SBATCH --job-name=FEniCSx_Array_Sims`**: Sets a collective name for the entire batch of 800 tasks.
  * **`#SBATCH --array=0-799`**: This is the crucial directive that creates the **Job Array**. It tells SLURM to generate 800 separate jobs (tasks), indexed from 0 up to and including 799.
  * **Resource Allocation (Per Task):**
      * `#SBATCH --nodes=1` and `#SBATCH --ntasks-per-node=1`: Each task runs as a single process (`mpirun`) on one node.
      * `#SBATCH --cpus-per-task=16`: Each of the 800 tasks is allocated **16 CPU cores** for its parallel execution (`mpirun`).
      * `#SBATCH --mem=64GB`: Each task is allocated **64GB of memory**.
  * **Output Management:**
      * `#SBATCH --output=./FEniCSx_3D_simulation_output/result_%A_%a.log`
      * `#SBATCH --error=./FEniCSx_3D_simulation_output/error_%A_%a.err`
      * The format specifiers `%A` (Job Array ID) and `%a` (Task ID, 0-799) ensure that each of the 800 tasks gets its own unique, non-overwriting output and error log file.

-----

### 🧮 Script Logic and Execution

1.  **Task Index Mapping:**

    ```bash
    NOK_VALUE=$SLURM_ARRAY_TASK_ID
    ```

    This line is key. The SLURM environment variable **`$SLURM_ARRAY_TASK_ID`** automatically contains the unique index of the current task (from 0 to 799). This index is stored in `NOK_VALUE` and is used to **differentiate the input** for each simulation (e.g., loading a specific parameter set for that index).

2.  **Environment Setup:**

    ```bash
    module purge
    module load fenicsx_complex/0.9.0
    ```

    The script ensures a clean environment and then loads the necessary **FEniCSx** module, which includes the required MPI implementation (likely Intel MPI in this case).

3.  **Parallel Execution:**

    ```bash
    mpirun -n $SLURM_CPUS_PER_TASK python ./<python_script.py> ${NOK_VALUE}
    ```

      * `mpirun -n $SLURM_CPUS_PER_TASK`: Initiates a parallel run using **16 processes** (the value of `$SLURM_CPUS_PER_TASK`).
      * `python ./<python_script.py>`: Executes the FEniCSx Python script.
      * `${NOK_VALUE}`: Passes the unique task index (0-799) as a **command-line argument** to the Python script. The Python script must be designed to read this argument and use it to select the correct parameter or input data for its specific run.

-----

### ➡️ Next Steps: Submitting the Job

To run all 800 simulations, you would submit this script using the standard SLURM command:

```bash
sbatch <script_name>.sh
```

SLURM will immediately create and queue 800 separate jobs based on the array definition.

#### Post-Analysis Dependency

The final section provides a powerful tip for workflow management:

```bash
sbatch --dependency=afterok:<Job_Array_ID> analyse_globale.sh
```

After submitting the array, you get a single **Job Array ID** (e.g., 123456). You can submit a separate global analysis script (`analyse_globale.sh`) and set a **dependency** on the array. This ensures the analysis script will **only start** after all 800 tasks in the array have completed successfully (`afterok`).

Would you like me to suggest how the Python script might use the `${NOK_VALUE}` argument?

### Sbatch file


```bash
#!/bin/bash

# ==============================================================================
# 1. SLURM DIRECTIVES (MODIFIED FOR JOB ARRAY 0-799)
# ==============================================================================
#SBATCH --job-name=FEniCSx_Array_Sims  # Main Job Array Name
#SBATCH --constraint=compute
#SBATCH --partition=i2m,i2m-resources
#
# --- JOB ARRAY PARAMETERS (800 tasks from 0 to 799) ---
#SBATCH --array=0-799
#
# --- RESOURCE PARAMETERS FOR EACH MPI TASK (Based on your successful test) ---
# Each task is an independent MPI run using 16 cores.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1             # One master process per node (mpirun)
#SBATCH --cpus-per-task=16              # The mpirun process will use 16 cores for MPI
#
# Time and Memory: Adjusted for 16 cores and 64GB (based on your test)
#SBATCH --time=0-01:00:00               # Time allocated PER TASK (adjust if necessary)
#SBATCH --mem=64GB                      # Memory allocated PER TASK (your test used ~46GB)
#
# Output Files: Slurm handles output redirection for you.
# %A = Job Array ID, %a = Task ID (0 to 799)
# The logs replace your previous './FEniCSx_3D_simulation_output/result${nok}.log' structure.
#SBATCH --output=./FEniCSx_3D_simulation_output/result_%A_%a.log
#SBATCH --error=./FEniCSx_3D_simulation_output/error_%A_%a.err
#
# --- OTHER PARAMETERS ---
#SBATCH --mail-type=END
#SBATCH --mail-user=<mail>@<domain>.fr # Replace with your email address
#SBATCH --chdir=/gpfs/home/<USERNAME>/your_working_dir # Set the working directory

# Exit immediately if a command exits with a non-zero status
set -e

# Slurm environment variable for the task index
NOK_VALUE=$SLURM_ARRAY_TASK_ID

# ==============================================================================
# 2. INFORMATION AND MODULE LOADING
# ==============================================================================
echo "==============================================================================" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID (Array ID):" $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID (NOK):" $NOK_VALUE
echo "Total processes (cores) requested:" $SLURM_CPUS_PER_TASK
echo "==============================================================================" 

echo "--- STARTING FENICSX SIMULATION #${NOK_VALUE} ---"

echo "Loading FEniCSx module (includes Python and Intel MPI)..."
# Ensure a clean environment before loading the FEniCSx module
module purge

# Load the installed FEniCSx environment module
module load fenicsx_complex/0.9.0

# ==============================================================================
# 3. EXECUTION 
# ==============================================================================

# Print job parameters for verification
echo "Starting execution on $(hostname) at $(date)"
echo "Number of cores requested: $SLURM_CPUS_PER_TASK"

# Execute the parallel FEniCSx script
# $SLURM_CPUS_PER_TASK is set to 16.
# The value $NOK_VALUE (from 0 to 799) is passed as an argument to the Python script.
echo "Running: mpirun -n $SLURM_CPUS_PER_TASK python ./<python_script.py> ${NOK_VALUE}"

mpirun -n $SLURM_CPUS_PER_TASK python ./<python_script.py> ${NOK_VALUE}

# ==============================================================================
# 4. POST-PROCESSING (OPTIONAL)
# ==============================================================================
echo "Job finished at $(date)"

echo "After a job finishes you can see the real use of CPUs and memory:"
echo "1. module load slurm/wrappers"
echo "2. seff <job-id>"
echo "3. module purge"

# ------------------------------------------------------------------------------
echo "To launch the global analysis job after all 800 simulations are complete:"
echo "Note the Job Array ID (e.g., 123456) after submitting this script with 'sbatch'."
echo "Submit a separate analysis script with a dependency:"
echo "sbatch --dependency=afterok:<Job_Array_ID> analyse_globale.sh"
# ------------------------------------------------------------------------------
```