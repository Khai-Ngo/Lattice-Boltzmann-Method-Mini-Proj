#!/bin/bash

#SBATCH --job-name=ompLBM
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --time=0:0:10
#SBATCH --mem=100M

# Load modules required for runtime
# no need to since I got it in .bashrc

cd $SLURM_SUBMIT_DIR

# Run the program
./ompLBM.exe 28 100 40 8500 25
