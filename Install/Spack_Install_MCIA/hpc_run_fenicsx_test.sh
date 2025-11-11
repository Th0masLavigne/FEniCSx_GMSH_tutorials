#!/bin/bash
#
# Slurm Script to test the FEniCSx/0.9.0 shared module and perform a custom pip install.
#
# ==============================================================================
# 1. SLURM DIRECTIVES
# ==============================================================================
#SBATCH --job-name=fenicsx_test_two_nodes
#SBATCH --constraint=compute
#SBATCH --nodes=2
#SBATCH --ntasks=4                      # Request 4 cores for MPI test
#SBATCH --tasks-per-node=2                # Request 4 cores per node
#SBATCH --time=0-00:35:00               # Short runtime
#SBATCH --mem=8GB                         # Use all memory on node
##SBATCH --mem-per-cpu=8g                # Memory per core
#SBATCH --output=slurm_test_%j.out
#SBATCH --error=slurm_test_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mail>@<domain>.fr
#SBATCH --chdir=/gpfs/home/<USERNAME>/your_working_dir

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

# Ensure the module environment is clean before loading
module purge

# Load the target module. This step sets $PATH, $PYTHONPATH, and loads Intel MPI.
module load fenicsx/0.9.0

# ==============================================================================
# 3. MPI EXECUTION TEST
# ==============================================================================

# Test the dolfinx import and version check using mpirun
# The 'mpirun' executable is provided by the Intel MPI module loaded by FEniCSx/0.9.0.
echo "Running parallel FEniCSx version check on $SLURM_NTASKS cores..."
echo "Elementary test..."

mpirun -n $SLURM_NTASKS python -c "from mpi4py import MPI; import dolfinx; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); size = comm.Get_size(); print(f'Hello from rank {rank} out of {size} processes, dolfinx v {dolfinx.__version__}')"

echo "From file test..."
mpirun -n $SLURM_NTASKS python Terzaghi.py

echo "--- TEST COMPLETE ---"