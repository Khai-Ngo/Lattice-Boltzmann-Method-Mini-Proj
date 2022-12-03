serial: 
	g++ -std=c++11 -o serialLBM.exe serialLBM.cpp -O3 -xHost
omp:
	g++ -std=c++11 -o ompLBM.exe openMP_LBM.cpp -O3 -xHost -fopenmp
all: serial omp
clean:
	rm -rf *.exe
