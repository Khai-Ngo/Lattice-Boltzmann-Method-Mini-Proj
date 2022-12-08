#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

    if (argc!=5){
        printf("Hey man, have 4 cml arguments\n");
        printf("<.exe><length><width><timeStepsNumber><saveFlag>\n");
        exit(1);
    }
    const int l = atoi(argv[1]);
    const int w = atoi(argv[2]);
    const int lat_size = l*w;
    const int time = atoi(argv[3]);
    const int saveFlag = atoi(argv[4]);
    float f[2][lat_size][9], rho[lat_size], xVel[lat_size], yVel[lat_size];

    printf("LBM of a %dx%d channel flow for %d_sec, west side inlet u0 = 0.2.\n. Serial C code.", l, w, time);
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
    int rowno, colno,e,west,n,s,ne,sw,nw,se;
    float dotProd, uProd, fiEQ_i;
    int z = 0; 
    clock_t begin = clock();
    for (int t=0;t<time;t++){
        for (int i=0; i<lat_size;i++){
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
            /* End of propagation phase */
            // Apply boundary conditions
            // left and right sidewalls
            for (int i = 0; i<lat_size;i+=w){
                // left bounce 5->6, 1->2, 8->7 
                f[1-z][i][6] = f[1-z][i][5];
                f[1-z][i][2] = f[1-z][i][1];
                f[1-z][i][7] = f[1-z][i][8];
                // right bounce 6->5, 2->1, 7->8
                f[1-z][i+w-1][5] = f[1-z][i+w-1][6];
                f[1-z][i+w-1][1] = f[1-z][i+w-1][2];
                f[1-z][i+w-1][8] = f[1-z][i+w-1][7];
            }
            // North src BC. 
            for (int i = 1; i < w-1; i++){
                density = 0;
                for (int j = 0; j < 9; j++){
                    density += f[1-z][i][j];
                }
                f[1-z][i][4] = f[1-z][i][2]+2*density*u0/3.0;
                f[1-z][i][6] = f[1-z][i][5] + 0.5*(f[1-z][i][1]-f[1-z][i][2])+density*u0/6.0;
                f[1-z][i][8] = f[1-z][i][7] - 0.5*(f[1-z][i][1]-f[1-z][i][2])+density*u0/6.0;
                }
            // South sink BC. 
            for (int i = lat_size - w; i< lat_size ; i++){
                f[1-z][i][3] = f[1-z][i-w][3];
                f[1-z][i][5] = f[1-z][i-w][5];
                f[1-z][i][7] = f[1-z][i-w][7];
            }             
            // update rho, xVel, yVel
            for (int i = 0; i < lat_size; i++){
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
    clock_t end = clock();
    double finalTime = (double)(end - begin) / CLOCKS_PER_SEC;
    // Write exports of rho and xVel map for visualization (courtesy of Python)
    printf("That took %f seconds\n", finalTime);
    printf("Received all data from workers. Now exporting...\n");
    if (saveFlag){
        char buffer[1024]; // export filename buffer
        snprintf(buffer, 1024, "rho_02u0_%dx%d_%dsec_serialC.txt", l,w,time);
        exportMap(buffer, l, w, &rho[0]);
        // snprintf(buffer, 1024, "xVel_02u0_%dx%d_%dsec_%dprocs.txt", l,w,time,numworkers);
        // exportMap(buffer,l,w,&xVel[0]);
        snprintf(buffer, 1024, "yVel_02u0_%dx%d_%dsec_serialC.txt", l,w,time);
        exportMap(buffer,l,w,&yVel[0]);
    }
    else
        printf("Actually user said no export\n");
    return 1;
}
