#!/bin/bash

#SBATCH --job-name=ompLBM
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=00:20:00
#SBATCH --mem=1000M
#SBATCH --array=100-1000:100

# Load modules required for runtime
# no need to since I got it in .bashrc

cd $SLURM_SUBMIT_DIR

# Run the program
./ompLBM.exe 28 $SLURM_ARRAY_TASK_ID 40 8500 0
