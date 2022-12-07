serial: 
	g++ -std=c++11 -o serialLBM.exe serialLBM.cpp -O3 -xHost
omp:
	g++ -std=c++11 -o ompLBM.exe openMP_LBM.cpp -O3 -xHost -fopenmp
mpi:
	mpicc -o mpi_LBM.exe mpi_LBM.c -O3 -xHost
serial_BC4:
	icpc -std=c++11 -o serialLBM_icc.exe serialLBM.cpp -O3 -xHost
omp_BC4:
	icpc -std=c+11 -o ompLBM_icc.exe open_LBM.cpp -O3 -xHost -qopenmp
mpi_BC4:
	mpiicc -o mpi_LBM_icc.exe mpi_LBM.c -O3 -xHost
all: serial omp mpi
all_BC4: serial_BC4 omp_BC4 mpi_BC4
clean:
	rm -rf *.exe
