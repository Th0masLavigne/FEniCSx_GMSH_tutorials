#!/bin/bash

# ==============================================================================
# 1. SLURM DIRECTIVES (MODIFY THESE FOR YOUR JOB)
# ==============================================================================
#SBATCH --job-name=fenicsx_convergence
#SBATCH --constraint=compute
#SBATCH --partition=i2m,i2m-resources
# --- JOB ARRAY PARAMETERS (4 tasks from 0 to 3) ---
#SBATCH --array=0-3 # for convergence analysis (i.e. use in mpirun) make it match --ntasks
#
#SBATCH --nodes=1 
#SBATCH --ntasks=4                      # Total number of MPI processes (4 total)
##SBATCH --ntasks-per-node=1              # One master process per node (mpirun)
##SBATCH --tasks-per-node=2              # Processes per node (2 per node on 2 nodes = 4 total)
#SBATCH --time=0-01:00:00               # Short runtime
#SBATCH --mem=16GB                       # Total memory reserved for the job
##SBATCH --mem-per-cpu=8g              # Alternative: Request memory per core
#SBATCH --output=/gpfs/home/tlavigne002/test_fenics/array/slurm_test_%j.out
#SBATCH --error=/gpfs/home/tlavigne002/test_fenics/array/slurm_test_%j.err
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<mail>@<domain>.fr # Replace with your email address
#SBATCH --chdir=/gpfs/home/tlavigne002/test_fenics/array

# Exit immediately if a command exits with a non-zero status
set -e

# Slurm environment variable for the task index (0-799)
SLURM_INDEX=$SLURM_ARRAY_TASK_ID

# Calculate the 1-based index (1-800)
# This is the value that will be passed to the Python script
NOK_VALUE=$((SLURM_INDEX + 1))

# useful informations to print
echo "==============================================================================" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID (Array ID):" $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_ID (0-799):" $SLURM_INDEX
echo "Calculated NOK Value (1-800):" $NOK_VALUE
echo "Total processes (cores) requested:" $SLURM_CPUS_PER_TASK
echo "==============================================================================" 

echo "--- STARTING FENICSX MODULE TEST #${NOK_VALUE} ---"

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
echo "Available MPI tasks: $SLURM_NTASKS"
echo "Total MPI tasks: $NOK_VALUE"

echo "Elementary test..."
mpirun -n ${NOK_VALUE} python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}');from petsc4py import PETSc; from dolfinx.fem.petsc import assemble_vector; print(PETSc.ScalarType);"


echo "--- STARTING FENICSX SIMULATION #${NOK_VALUE} ---"
# Execute the parallel Python script using mpirun
# mpirun uses the OpenMPI instance loaded by the FEniCSx module
echo "Running: mpirun -n $NOK_VALUE python Hyper_elastic_contact_conditionnel_body_force.py "

# IMPORTANT: Ensure <filename>.py is in the directory specified by --chdir, add the nok_value for your export names for instance
mpirun -n ${NOK_VALUE} python Hyper_elastic_contact_conditionnel_body_force.py ${NOK_VALUE}

# ==============================================================================
# 4. POST-PROCESSING (OPTIONAL)
# ==============================================================================
echo "Job finished at $(date)"

echo "After a job finishes you can see the real use of CPUs and memory"
echo "1. module load slurm/wrappers"
echo "2. seff <job-id>"
echo "3. module purge"