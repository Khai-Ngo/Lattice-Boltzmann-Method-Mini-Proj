# Lattice Boltzmann Method Mini Project
Advanced Computational Physics mini project at the University of Bristol. The aim is to investigate how the performance of a Lattice Boltzmann code scales with degreee of parallelisation and problem size. Fluid flow through a narrow channel is simulated and timed as the means of investigation.  

## Compiliation instruction: 

Serial LBM code: 

make serial -> builds the OOP (i.e. slower) version of the serial code using gcc
make serialC_BC4 -> builds straight C version of the serial code using icc (for use on BC4)
make serial_BC4 -> builds the OOP version using icc compiler for BC4 runs

OpenMP code:

make omp -> builds the OOP code parallelised using openMP
make omp_BC4 -> same as above, but using icc for BC4 runs

MPI code: 

make mpi -> builds the C code parallelised using OpenMPI
make mpi_BC4 -> same as above, but using icc for BC4 runs

## About command line arguments:

The serial code takes in 4 arguments
<lattice_length><lattice_width><number_of_time_steps><save_flag>
(save_flag = 1 saves the final density and horizontal velocity map in 2 text files. Set this to 0 if printing is not desired.)

The openMP code takes in 5 arguments
<max_number_of_threads><lattice_length><lattice_width><number_of_time_steps><save_flag>
(save_flag = n saves the final density and horizontal velocity map for the run that used n threads in 2 text files. Set this to 0 if printing is not desired.)

The MPI code takes in 4 arguments exactly the same way as the serial code

For example: 

./ompLBM.exe 10 1000 40 8500 5

tells the program to compute the channel flow for a 1000x40 grid, for 8500 steps, for 10 times with the number of threads set to 1, 2, 3, ... , with execution times taken for each iteration. The density and horizontal velocity maps are written into file for the run with 5 threads. 



