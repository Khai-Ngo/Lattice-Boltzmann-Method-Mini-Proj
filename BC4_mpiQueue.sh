#!/bin/bash

#SBATCH --job-name=mpiLBM
#SBATCH --partition=teach_cpu
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=1000M
#SBATCH --array=3-28:1

# Load modules required for runtime
# no need to since I got it in .bashrc

cd $SLURM_SUBMIT_DIR

# Run the program
mpiexec -n $SLURM_ARRAY_TASK_ID ./mpi_LBM_icc.exe 100 40 8500 0
