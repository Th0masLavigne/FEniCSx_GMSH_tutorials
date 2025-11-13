#!/bin/bash

# ==============================================================================
# 1. SLURM DIRECTIVES (MODIFY THESE FOR YOUR JOB)
# ==============================================================================
#SBATCH --job-name=fenicsx_test_two_nodes
#SBATCH --constraint=compute
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

echo "After a job finishes you can see the real use of CPUs and memory"
echo "1. module load slurm/wrappers"
echo "2. seff <job-id>"
echo "3. module purge"