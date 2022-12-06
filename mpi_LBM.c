#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>

// This code solves flow along channel on a transposed matrix (i.e. a vertical channel flowing from north to south) so parallelization is convenient (maybe better performance?)
// so you have left and right (otherwise it's above and below)
// python shall do the magic of transposing the matrix and plot it as flow from left to right again 
// also obvious that it's based on heat equation code 

void initLattice(int lat_size, float *f){
    // initialize the lattice in master task. Note f[1] is not touched since there's no need for it. 
    // master shall later received results from each worker and paste it into this initialized lattice
    for (int i = 0; i <lat_size ;i++){
            *(f+i*9) = 1; 
            for (int k =1 ; k< 9; k++){
                *(f+i*9+k) = 0;
            } 
    }       
}
void updateRUV(int lat_size, int* ex, int* ey, float *rho, float *xVel, float *yVel, float *f){
    // probably to be run in the master task only
    float vx, vy, density;
    for (int i = 0; i < lat_size; i++)
    {
        vx = 0;
        vy = 0;
        density =0; 
        for (int k =0; k<9; k++){
            vx+= ex[k]*(*(f+i*9+k));
            vy+= ey[k]*(*(f+i*9+k));
            density+= (*(f+i*9+k));
        }
        vx /= density;
        vy /= density;
        *(rho+i) = density;
        *(xVel+i) = vx;
        *(yVel+i) = vy;
    }
}
void exportMap(char* fname, int l, int w, float* map){
    // visualizer programm shall need to transpose this
    int lat_size = l*w;
    FILE *f;
    f = fopen(fname, "w");
    for (int i = 0; i<lat_size; i++){
            char* endchar = ((i%w)==(w-1))? "\n":"\t";
            fprintf(f, "%8.1f%s", *(map+i), endchar);

    }
    fclose(f);
    printf("Successfully wrote data to %s", fname);
}

int main(int argc, char* argv[]){
    void initLattice(), ruv(), exportMap();
    // MPI related variables
    int taskid, numworkers, numtasks, 
    averow, rows, offset, extra,
    dest, src, 
    left, right,
    msgtype, 
    errCode, start,end;
    MPI_Status status;

    const unsigned MASTER=0;
    const unsigned NONE=0; // no neighbour 
    const unsigned MAXWORKER=27; // work divided at most 27 BC4 cores. Leave 1 core to be the boss
    const unsigned MINWORKER=2; // work divided at least between 2 cores. Third core being the boss.
    const unsigned BEGIN=1; // msg tag
    const unsigned DONE=4; // msg tag
    const unsigned LTAG=2;//msg tag
    const unsigned RTAG=3; //msg tag
    // LBM related variables
    // North side source (should be West when matrix is transposed)
    const float u0=0.2;
    const float alpha=0.02;
    const float tau = 3.0*alpha+0.5; 
    // the weights and ei vectors in D2Q9 scheme
    const float weights[9] = {4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};
    unsigned ex[9] = {0,1,-1,0,0,1,-1,-1,1}; // hope to god is not overwritten LOL
    unsigned ey[9] = {0,0,0,1,-1,1,-1,1,-1};

    // find out how many tasks are running, and your current taskid
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers = numtasks-1;
    // MASTER TASK DO THE BELOW ONLY:
    if (taskid == MASTER){
        /* only master checks for cml argument*/
    if (argc!=5){
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
    const int lat_size = l*w;
    const int time = atoi(argv[3]);
    const int saveFlag = atoi(argv[4]);
    float f[lat_size][9], rho[lat_size], xVel[lat_size], yVel[lat_size];
    // and only the boss prints stuff
    printf("%dx%d_%d_sec\n", l, w, time);
    initLattice(lat_size, &f[0][0]); // note f stores the address of f[0][0][0] (i.e. f_0 at (0, 0) of the current lattice. f[1][0][0] is f_0 at (0,0) for the duplicate lattice)
    updateRUV(lat_size, ex, ey, &rho[0], &xVel[0], &yVel[0], &f[0][0]);
    //exportMap("init_rho.txt", l, w, &rho[0]); //sanity check
    //exportMap("init_xVel.txt", l, w, &xVel[0]);  //sanity check
    
    // start timing from when you start sending work (all the f_i of every lattice point, and the array partitions)
    int startTime, finalTime;
    startTime = MPI_Wtime();
    // Distribute work
    averow = l/numworkers; // number of row per worker
    extra = l%numworkers; // extra rows. A number between 0 and numworkers-1
    offset=0; 
    for (unsigned i =1 ; i<numworkers+1;i++){
        dest = i;
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
        MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&f[offset][0], rows*w, MPI_FLOAT, dest, BEGIN,MPI_COMM_WORLD);
        printf("Sent to task %d: rows= %d offset= %d ",dest,rows,offset);
        printf("left= %d right= %d\n",left,right);
        offset = offset + rows*w;
    }
    // now wait for work done by loyal workers
    for (unsigned i = 1; i< numworkers+1;i++){
        src = i;
        msgtype = DONE;
        MPI_Recv(&offset, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&f[offset][0], rows*w, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
    }
    // stop timing when workers have received, done its share of work, and sent back all f_i data to master
    finalTime = MPI_Wtime() - startTime;

    // Write exports of rho and xVel map for visualization (courtesy of Python)
    updateRUV(lat_size, ex, ey, &rho[0], &xVel[0], &yVel[0], &f[0][0]); // I don't include this in code timing because I didn't either for OpenMP code. Problem is consider solved when all f_i are solved after t timesteps.
    printf("Received all data from workers. Now exporting...\n");
    char buffer[1024]; // export filename buffer
    snprintf(buffer, 1024, "rho_02u0_%dx%d_%dsec_%dprocs.txt", l,w,time,numworkers);
    exportMap(buffer, l, w, &rho[0]);
    snprintf(buffer, 1024, "xVel_02u0_%dx%d_%dsec_%dprocs.txt", l,w,time,numworkers);
    exportMap(buffer,l,w,&xVel[0]);
    }
    // WORKERS DO THE BELOW
    if (taskid!=MASTER){
        // worker creates the same grid lattice, but they only work on their given portion
          
    }
    return 1;
}
