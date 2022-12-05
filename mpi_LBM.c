#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

void collide(double *f);
void exportRho(char* fname, double* rho[]);
void exportxVel(char* fname, double* xVel[]);

int main(int argc, char* argv[]){
    // MPI related variables
    int taskid, numworkers, numtasks, 
    averow, rows, offset, extra,
    dest, src, 
    left, right,
    msgtype, 
    errCode, start,end;
    MPI_Status status;

    const unsigned MASTER=0;
    const unsigned NONE=0;
    const unsigned MAXWORKER=27; // work divided at most 27 BC4 cores. Leave 1 core to be the boss
    const unsigned MINWORKER=3; // work divided at least between 2 cores. Third core being the boss.
    // LBM related variables
    // West side source
    const double u0=0.2;
    const double alpha=0.02;
    const double tau = 3.0*alpha+0.5; 
    // the weights and ei vectors in D2Q9 scheme
    const double weights[9] = {4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};
    const unsigned ex[9] = {0,1,-1,0,0,1,-1,-1,1};
    const unsigned ey[9] = {0,0,0,1,-1,1,-1,1,-1};

    // find out how many tasks are running, and your current taskid
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers = numtasks-1;
    // MASTER TASK DO THE BELOW ONLY:
    if (taskid == MASTER){
        /* only master checks for cml argument*/
    if (argc!=6){
        printf("Hey man, have 4 cml arguments\n");
        printf("<.exe><length><width><timeStepsNumber><saveFlag>\n");
        MPI_Abort(MPI_COMM_WORLD, errCode);
        exit(1);
    }
    if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)){
        printf("Hey man, too many or too few workers!\n");
        printf("I quit\n");
        MPI_Abort(MPI_COMM_WORLD, errCode);
        exit(1);
    }
    // boss defines the grid and stuff
    const int l = atoi(argv[1]);
    const int w = atoi(argv[2]);
    const int time = atoi(argv[3]);
    const int saveFlag = atoi(argv[4]);
    double f_cur[l][w][9],  f_dup[l][w][9], rho[l][w], xVel[l][w];
    // and only the boss prints stuff
    printf("%dx%d_%d_sec\n", l, w, time);
    // Distribute work
    averow = l/numworkers; // number of row per worker
    extra = l%numworkers; // extra rows. A number between 0 and numworkers-1
    offset=0; 
    for (int i =1 ; i<numworkers+1;i++){
        // each task take on 1 extra row, until you've got no extra row (worst case everyone took on an extra row except for 1 task)
        rows=(i<=extra)? averow+1:averow; 
        // Tell workers who their neighbours are
        if (i==1)
            left = NONE;
        else 
            left = i-1;
        if (i == numworkers)
            right = NONE;
        else right = i+1;
        // Send startup info to each worker
        dest = i;
        

    }
    }
    // WORKERS DO THE BELOW
    if (taskid!=MASTER){
        // worker creates the same grid lattice, but they only work on their given portion

    }
    return 1;
}
