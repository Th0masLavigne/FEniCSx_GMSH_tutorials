#!/bin/bash

#SBATCH --job-name=ID_63_refine_para
#SBATCH --partition=CPU_Compute
#SBATCH --ntasks=8
#SBATCH --nodes=1
##SBATCH --tasks-per-node=1
##SBATCH --mem-per-cpu=32GB
#SBATCH --mem=0
#
## Suggested batch arguments
##SBATCH --mail-type=ALL
##SBATCH --mail-user=thomas.lavigne@ensam.eu
#
## Logging arguments (IMPORTANT)
#SBATCH --output=slurm_%x-%j.out
#SBATCH --error=slurm_%x-%j.err


. /modules/spack/v0.23/share/spack/setup-env.sh
spack env activate tlavigne

# srun python main_finger.py
mpirun -n $SLURM_NTASKS python main_finger.py

