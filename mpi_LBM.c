#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>

// This code solves flow along channel on a transposed matrix (i.e. a vertical channel flowing from north to south) so parallelization is convenient (maybe better performance?)
// left and right should be understood as aliases for above and below in this code. 
// python shall do the magic of transposing the matrix and plot it as flow from left to right again 
// also obvious that it's based on heat equation code 


void exportMap(char* fname, int l, int w, float* map){
    // visualizer programm shall need to transpose this
    int lat_size = l*w;
    FILE *f;
    f = fopen(fname, "w");
    for (int i = 0; i<lat_size; i++){
            char* endchar = ((i%w)==(w-1))? "\n":"\t";
            fprintf(f, "%f%s", map[i], endchar);
    }
    fclose(f);
    printf("Successfully wrote data to %s\n", fname);
}

int main(int argc, char* argv[]){
    void exportMap();
    // MPI related variables
    int taskid, numworkers, numtasks, 
    aveBlock, chunks, offset, extra, averow,
    dest, src, 
    left, right,
    msgtype, 
    errCode, start,end;
    MPI_Status status;

    const unsigned MASTER=0;
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
    const float tau = 3.0*alpha+0.5; // 0.56
    // temp vars
    float vx, vy, density;
    
    // the weights and ei vectors in D2Q9 scheme
    float weights[9] = {4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};
    int ex[9] = {0,1,-1,0,0,1,-1,-1,1}; // hope to god is not overwritten since if passed to function you can't have const unsigned
    int ey[9] = {0,0,0,1,-1,1,-1,1,-1};

    // find out how many tasks are running, and your current taskid
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    numworkers = numtasks-1;
    // brave, but at least known to all tasks from the start
    const int l = atoi(argv[1]);
    const int w = atoi(argv[2]);
    const int lat_size = l*w;
    const int time = atoi(argv[3]);
    const int saveFlag = atoi(argv[4]);
    float f[2][lat_size][9], rho[lat_size], xVel[lat_size], yVel[lat_size];

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
    // only master prints stuff
    printf("LBM of a %dx%d channel flow for %d_sec, west side inlet u0 = 0.2. %d worker tasks\n", l, w, time, numworkers);
    // lattice intialization in master
    for (int k =0; k < 2; k++){
        for (int i=0;i<lat_size;i++){
            f[k][i][0]=1;
            for (int j=1;j<9;j++){
                f[k][i][j]=0;
            }
        }
    }
    // first calc of rho, xVel, yVel maps right after init
    for (int i = 0; i < lat_size; i++){
        vx = 0;
        vy = 0;
        density = 0;
        for (int j = 0; j < 9; j++){
            vx+=ex[j]*f[0][i][j];
            vy += ey[j]*f[0][i][j];
            density += f[0][i][j];
        }
        rho[i] = density;
        xVel[i] = vx/density;
        yVel[i] = vy/density;
    }
    
    // start timing from when you start sending work (all the f_i of every lattice point, and the array partitions)
    // cheat a bit by not including lattice initialization
    double startTime, finalTime;
    startTime = MPI_Wtime();
    // Now you distribute work
    averow = l/numworkers; // number of row per worker
    extra = l%numworkers; // extra rows. A number between 0 and numworkers-1
    offset=0; 
    for (unsigned i =1 ; i<numworkers+1;i++){
        dest = i;
        // each task take on 1 extra row, until you've got no extra row (worst case everyone took on an extra row except for 1 task)
        chunks= (i<=extra)? (averow+1)*w:averow*w; 
        // Tell workers who their neighbours are
        if (i==1)
            left = numworkers; // left flow loops all the way to the back 
        else 
            left = i-1;
        if (i == numworkers)
            right = 1; // right flow loops all the way to the front
        else right = i+1;
        // Send startup info to each worker
        MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&chunks, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&numworkers, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&f[0][offset][0], chunks*9, MPI_FLOAT, dest, BEGIN,MPI_COMM_WORLD);
        MPI_Send(&rho[offset], chunks, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&xVel[offset], chunks, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);
        MPI_Send(&yVel[offset], chunks, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);
        printf("Sent to task %d: chunks= %d offset= %d ",dest,chunks,offset);
        printf("left= %d right= %d\n",left,right);
        offset += chunks;
    }
    // now wait for work done by loyal workers
    for (unsigned i = 1; i< numworkers+1;i++){
        src = i;
        msgtype = DONE;
        MPI_Recv(&offset, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&chunks, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&f[0][offset][0], chunks*9, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rho[offset], chunks, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&xVel[offset], chunks, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&yVel[offset], chunks, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
    }
    // stop timing when workers have received, done its share of work, and sent back all f_i data to master
    finalTime = MPI_Wtime() - startTime;

    // Write exports of rho and xVel map for visualization (courtesy of Python)
    printf("That took %f seconds\n", finalTime);
    printf("Received all data from workers. Now exporting...\n");
    if (saveFlag){
        char buffer[1024]; // export filename buffer
        snprintf(buffer, 1024, "rho_02u0_%dx%d_%dsec_%dprocs.txt", l,w,time,numworkers);
        exportMap(buffer, l, w, &rho[0]);
        // snprintf(buffer, 1024, "xVel_02u0_%dx%d_%dsec_%dprocs.txt", l,w,time,numworkers);
        // exportMap(buffer,l,w,&xVel[0]);
        snprintf(buffer, 1024, "yVel_02u0_%dx%d_%dsec_%dprocs.txt", l,w,time,numworkers);
        exportMap(buffer,l,w,&yVel[0]);
    }
    else
        printf("Actually user said no export\n");
    MPI_Finalize();
    // end of master code
    }
    // WORKERS DO THE BELOW
    if (taskid!=MASTER){
        // worker creates the same grid lattice, but they only work on their given portion
       // Initialize lattice just like master first
       for (int k =0; k < 2; k++){
        for (int i=0;i<lat_size;i++){
            f[k][i][0]=1;
            for (int j=1;j<9;j++){
                f[k][i][j]=0;
                }
            }
        }
        // receive its share of the lattice
        src = MASTER;
        msgtype = BEGIN;
        MPI_Recv(&offset, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&chunks, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&left, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&right, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&numworkers, 1, MPI_INT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&f[0][offset][0], chunks*9, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rho[offset], chunks, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&xVel[offset], chunks, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&yVel[offset], chunks, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
        // process its share of the lattice 
        start = offset;
        end = offset+chunks;
        int rowno, colno,e,west,n,s,ne,sw,nw,se;
        float dotProd, uProd, fiEQ_i;
        int z = 0; 
        for (int t=0;t<time;t++){
            for (int i=start; i<end;i++){
            // collision
            for (int j = 0;j<9;j++){
                dotProd = ex[j]*xVel[i]+ey[j]*yVel[i];
                uProd = xVel[i]*xVel[i] + yVel[i]*yVel[i];
                fiEQ_i = weights[j]*rho[i]*(1+3*dotProd+4.5*dotProd*dotProd-1.5*uProd);
                f[z][i][j] -= (f[z][i][j]-fiEQ_i)/tau;
            }
            /* Start of propagation phase */
            rowno = i/w;
            colno = i%w;
            e = ((colno != (w-1)))? i+1 : i+1-w;
            west = ((colno != 0))? i-1:i-1+w;
            n =  (rowno != 0) ? i-w : i - w + lat_size; 
            s = (rowno != (l-1))? i+w : i + w - lat_size;
            ne = ((colno != (w-1)))? i-w+1 : i-2*w+1;
            if (ne<0) ne+= lat_size;
            sw = ((colno != 0))? i+w-1:i+2*w-1;
            if (sw >=lat_size) sw-=lat_size;
            nw = ((colno != 0))? i-w-1:i-1;
            if(nw<0) nw+= lat_size;
            se = ((colno != (w-1)))? i+w+1 : i+1;
            if (se >= lat_size) se-=lat_size;
            f[1-z][i][0] = f[z][i][0]; // e0
            f[1-z][e][1] = f[z][i][1]; // e1 E
            f[1-z][west][2] = f[z][i][2]; // e2 W 
            f[1-z][n][3] = f[z][i][3]; // e3 N
            f[1-z][s][4] = f[z][i][4]; // e4 S
            f[1-z][ne][5] = f[z][i][5]; // e5 NE
            f[1-z][sw][6] = f[z][i][6]; // e6 SW
            f[1-z][nw][7] = f[z][i][7]; // e7 NW
            f[1-z][se][8] = f[z][i][8]; // e8 SE
            }
            // communicate with neighboring blocks
            // send row before 1st row of this block (that this block propagated into) to above, or last row to the last block if left == numworkers
            if (left == numworkers){
                MPI_Send(&f[1-z][lat_size-w][0], w*9, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
            }
            else{
                MPI_Send(&f[1-z][start - w][0], w*9, MPI_FLOAT, left, RTAG, MPI_COMM_WORLD);
            }

            // send row after last row of this block (that this block propagated into) to below, or first row to the first block if right == 1
            if (right == 1){
                MPI_Send(&f[1-z][0][0], w*9, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);

            }
            else{
                MPI_Send(&f[1-z][end][0], w*9, MPI_FLOAT, right, LTAG, MPI_COMM_WORLD);
            }
            
            // receive data propagated from below up here, or from the first box onto the last row
            src=right;
            msgtype=RTAG;
            MPI_Recv(&f[z][end-w][0], w*9, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
            // then update e3, e5, e7
            for (int i=end-w; i<end;i++){
                f[1-z][i][3] = f[z][i][3];
                f[1-z][i][5] = f[z][i][5];
                f[1-z][i][7] = f[z][i][7];
            }
            // receive data propagated from above down here, or from the last box onto the first row
            src = left;
            msgtype = LTAG;
            MPI_Recv(&f[z][start][0],w*9, MPI_FLOAT, src, msgtype, MPI_COMM_WORLD, &status);
            // then update e4, e6, e8
            for (int i= start; i<start+w;i++){
                f[1-z][i][4] = f[z][i][4];
                f[1-z][i][6] = f[z][i][6];
                f[1-z][i][8] = f[z][i][8];
            }
            /* End of propagation phase */
            // Apply boundary conditions
            // left and right sidewalls
            for (int i = start; i<end;i+=w){
                // left bounce 5->6, 1->2, 8->7 
                f[1-z][i][6] = f[1-z][i][5];
                f[1-z][i][2] = f[1-z][i][1];
                f[1-z][i][7] = f[1-z][i][8];
                // right bounce 6->5, 2->1, 7->8
                f[1-z][i+w-1][5] = f[1-z][i+w-1][6];
                f[1-z][i+w-1][1] = f[1-z][i+w-1][2];
                f[1-z][i+w-1][8] = f[1-z][i+w-1][7];
            }
            // North src BC. This will be the bottle neck since only task 1 does it. Requires better load balancing.
            if (start == 0){
                for (int i = 1; i < w-1; i++){
                    density = 0;
                    for (int j = 0; j < 9; j++){
                        density += f[1-z][i][j];
                    }
                     f[1-z][i][4] = f[1-z][i][2]+2*density*u0/3.0;
                    f[1-z][i][6] = f[1-z][i][5] + 0.5*(f[1-z][i][1]-f[1-z][i][2])+density*u0/6.0;
                    f[1-z][i][8] = f[1-z][i][7] - 0.5*(f[1-z][i][1]-f[1-z][i][2])+density*u0/6.0;
                    }
            }
            // South sink BC. This will be the bottle neck since only the last task does it. Requires better load balancing.
            if (end == lat_size){
                for (int i = end - w; i< end ; i++){
                    f[1-z][i][3] = f[1-z][i-w][3];
                    f[1-z][i][5] = f[1-z][i-w][5];
                    f[1-z][i][7] = f[1-z][i-w][7];
                }             
            }
            // update rho, xVel, yVel
            for (int i = start; i < end; i++){
                    vx = 0;
                    vy = 0;
                    density = 0;
                    for (int j = 0; j < 9; j++){
                        vx += ex[j]*f[1-z][i][j];
                        vy += ey[j]*f[1-z][i][j];
                        density += f[1-z][i][j];
                    }
                    rho[i] = density;
                    xVel[i] = vx/density;
                    yVel[i] = vy/density;
                    }
            z = 1-z; // curent set is the new set
        }
        // finally send work back to master
        MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&chunks, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&f[z][offset][0], chunks*9, MPI_FLOAT, MASTER, DONE,MPI_COMM_WORLD);
        MPI_Send(&rho[offset], chunks, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&xVel[offset], chunks, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Send(&yVel[offset], chunks, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
        MPI_Finalize();
    }
    return 1;
}
