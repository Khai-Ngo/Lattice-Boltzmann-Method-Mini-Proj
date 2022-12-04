#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char* argv[]){
    
    const double u0=0.2;
    const double alpha=0.02;
    const double tau = 3.0*alpha+0.5; 

    const double weights[9] = {4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};
    const unsigned ex[9] = {0,1,-1,0,0,1,-1,-1,1};
    const unsigned ey[9] = {0,0,0,1,-1,1,-1,1,-1};

    if (argc!=6){
        printf("Hey man, 5 cml arguments\n");
        printf("<.exe><maxProcNo><length><width><timeStepsNumber><saveFlag>\n");
        return 0;
    }

    const int maxProcNo = atoi(argv[1]);
    const int l = atoi(argv[2]);
    const int w = atoi(argv[3]);
    const int time = atoi(argv[4]);
    const int saveFlag = atoi(argv[5]);
    double results[maxProcNo];
    double f_cur[l][w][9], f_new[l][w][9], rho[l][w], xVel[l][w], yVel[l][w];
    printf("%dx%d_%d_sec\n", l, w, time); 
    return 1;
}
