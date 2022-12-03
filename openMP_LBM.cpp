#include <vector>
#include <array>
#include <iostream>
#include <utility>
#include <fstream>
#include <cmath>
#include <string>
#include <omp.h>

namespace LBM{
class Cell{
    private:
        const std::array<double, 9> weights {{4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36}};
        // like the diagram in the lecture notes: static, E, W, N, S, NE, SW, NW, SE
        const std::array<int, 9>ex {{0,1,-1,0,0,1,-1,-1,1}};
        const std::array<int, 9>ey {{0,0,0,1,-1,1,-1,1,-1}};
        std::array<double, 9>f_i;
        std::array<double, 9>f_dup; // for updating. Avoid overwriting data before they get propagated
        double rho;
        double xVel;
        double yVel;
        void ruv(){
            double vx=0;
            double vy=0;
            double density=0;
            for (int i=0;i<f_i.size();i++){
                vx+=ex[i]*f_i[i];
                vy+=ey[i]*f_i[i];
                density+=f_i[i]; 
            }
            vx /= density;
            vy /= density;
            this->rho = density;
            this->xVel = vx;
            this->yVel= vy;
        }
    public:
        Cell(){
             // concurrent f_i and its duplicate. Default to rho=1, no velocity
            this->f_i = {1, 0, 0, 0, 0, 0, 0, 0, 0};
            this->f_dup = {1, 0, 0, 0, 0, 0, 0, 0, 0};
            ruv();
        }
        // copy constructor
        Cell (const Cell& inputCell){
            this->f_i = inputCell.f_i;
            this->f_dup = inputCell.f_dup;
            ruv();
        }
        void update_fi(int dir, double input){
            this->f_dup[dir] = input;
        }
        double get_fi(int dir){
            return this->f_i[dir];
        }
        void sync(){
            this->f_i = this->f_dup;
            ruv();
        }
        // Avoid looping through the array many times but repeat code
        void collision(double t){
            double dotProd, uProd, fiEQ_i;
            for (int i = 0; i<f_i.size();i++){  
                dotProd = ex[i]*xVel+ey[i]*yVel;
                uProd = xVel*xVel+yVel*yVel;
                fiEQ_i= weights[i]*rho*(1+3*dotProd+4.5*dotProd*dotProd-1.5*uProd);
                f_i[i] -= (f_i[i]-fiEQ_i)/t;
            }
        }
        double density(){
            return rho;
        }
        double vx(){
            return xVel;
        }
        double vy(){
            return yVel;
        }
        void cellReset(){
            this->f_i = {1, 0, 0, 0, 0, 0, 0, 0, 0};
            this->f_dup = {1, 0, 0, 0, 0, 0, 0, 0, 0};
            ruv();
        }
        // for debugging purposes
        std::array<double, 9> distFuctions(){
            return std::move(f_i);
        }
        std::array<double, 9> equilibrium(){
            std::array<double, 9> fiEQ;
            for (int i = 0; i<f_i.size();i++){  
                double dotProd = ex[i]*xVel+ey[i]*yVel;
                double uProd = xVel*xVel+yVel*yVel;
                fiEQ[i]= weights[i]*rho*(1+3*dotProd+4.5*dotProd*dotProd-1.5*uProd);
            }
            return std::move(fiEQ); //  to avoid memory run-away
        }
};
class LatticeBoltzmann{
    private:
        int l;
        int w;
        int lat_size;
        double tau;
        std::vector<Cell> cellLattice;
        int xy_to_index(int x, int y, int l){
            return x*l+y;
        }
        // private helper functions
        // note: three are 3 cases: top, bottom, and middle
        bool isTop(int index){
            return (index/l==0); // index/l is row number in case of 2D indexing
        }
        bool isBottom(int index){
            return(index/l==(w-1));
        }
        // also 3 cases: left, right, and middle
        bool isLeft(int index){
            return (index%l==0); // index%l is column number in case of 2D indexing
        }
        bool isRight(int index){
            return (index%l==(l-1));
        }
        inline void syncAll(){
            #pragma omp parallel for
            for (int j = 0; j<lat_size;j++){
                cellLattice[j].sync();
            } 
        }
    public:
        LatticeBoltzmann(){};
        LatticeBoltzmann(int l, int w, double tau){
            this->l = l;
            this->w = w;
            this->tau = tau;
            this->lat_size = l*w;
            // initialize all lattice cells to rho=1 everywhere. No flow. 
            Cell *initCell = new Cell();         
            for (int i=0; i<lat_size; i++){        
                this->cellLattice.push_back(*initCell);
            }
            delete initCell;      
        }
        double simulate(double u0, int time, int nThreads){
            omp_set_num_threads(nThreads);
            int e,w,n,s,ne,sw,nw,se;
            double rho, f1update, f5update, f8update;
            double initial, final;

            initial=omp_get_wtime();

            for (int t=0; t<time;t++){
                #pragma omp parallel private(e,w,n,s,ne,sw,nw,se)
                {
                #pragma omp for
                    for (int i=0;i<lat_size;i++){
                    // collision step                 
                    this->cellLattice[i].collision(tau);
                    // advection step, with periodic looping around
                    e = (!isRight(i))? i+1 : i+1-l;
                    w = (!isLeft(i))? i-1 : i-1+l;
                    n = (!isTop(i))? i-l : i - l+ lat_size;
                    s = (!isBottom(i))? i+l : i + l - lat_size;
                    ne = (!isRight(i))? i-l+1 : i-2*l+1;
                    if (ne<0) ne+=lat_size;
                    sw = (!isLeft(i))? i+l-1 : i+2*l-1;
                    if (sw>=lat_size) sw-=lat_size;
                    nw = (!isLeft(i))? i-l-1 : i -1;
                    if (nw <0) nw+=lat_size;
                    se = (!isRight(i))? i+l+1 : i+1;
                    if (se >= lat_size) se-=lat_size;
                
                    this->cellLattice[e].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
                    this->cellLattice[w].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
                    this->cellLattice[n].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
                    this->cellLattice[s].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
                    this->cellLattice[ne].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE
                    this->cellLattice[sw].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW
                    this->cellLattice[nw].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW
                    this->cellLattice[se].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE
                    
                    }
                }  
                // sync all propagated fluxes and update all rho and velocities
                syncAll();
                // Apply boundary condition
                #pragma omp parallel for 
                for (int i=0; i<l;i++){
                    // top bounce
                    this->cellLattice[i].update_fi(3, cellLattice[i].get_fi(4));
                    this->cellLattice[i].update_fi(5, cellLattice[i].get_fi(6));
                    this->cellLattice[i].update_fi(7, cellLattice[i].get_fi(8));
                    this->cellLattice[i].sync();
                    // bottom bounce
                    this->cellLattice[i+lat_size-l].update_fi(4, cellLattice[i+lat_size-l].get_fi(3));
                    this->cellLattice[i+lat_size-l].update_fi(6, cellLattice[i+lat_size-l].get_fi(5));
                    this->cellLattice[i+lat_size-l].update_fi(8, cellLattice[i+lat_size-l].get_fi(7));
                    this->cellLattice[i+lat_size-l].sync();
                }                    
                // inlet fixed x-velocity
                // this step actually also initializes the west side inlet flow 
                #pragma omp parallel private(rho, f1update, f5update, f8update) 
                #pragma omp for 
                for (int j=l;j<=lat_size-2*l;j+=l){
                    // west side sink
                    rho = cellLattice[j].density();
                    f1update= cellLattice[j].get_fi(2)+2*rho*u0/3.0;
                    f5update= cellLattice[j].get_fi(6)-0.5*(cellLattice[j].get_fi(3)-cellLattice[j].get_fi(4))+rho*u0/6.0;
                    f8update = cellLattice[j].get_fi(7)+0.5*(cellLattice[j].get_fi(3)-cellLattice[j].get_fi(4))+rho*u0/6.0;
                    this->cellLattice[j].update_fi(1, f1update);
                    this->cellLattice[j].update_fi(5, f5update);
                    this->cellLattice[j].update_fi(8, f8update);
                    this->cellLattice[j].sync();
                }
                #pragma omp parallel for
                for (int j=l-1;j<lat_size;j+=l){
                    this->cellLattice[j].update_fi(2, cellLattice[j-1].get_fi(2));
                    this->cellLattice[j].update_fi(7, cellLattice[j-1].get_fi(7));
                    this->cellLattice[j].update_fi(6, cellLattice[j-1].get_fi(6));
                    this->cellLattice[j].sync();
                }
            }
            final = omp_get_wtime();
            return final-initial;
        }
        void exportRho(std::string& fname){
            std::ofstream f (fname);
            if (!f){
                std::cout<<"Can't open "<<fname<<" for writing!\n";
            }
            for (int i=0; i<lat_size;i++){
                double rho = cellLattice[i].density();
                std::string endChar = (i%l==(l-1))? "\n":"\t";
                f<<rho<<endChar;
            }
            f.close();
        }
        void exportxVel(std::string& fname){
            std::ofstream f (fname);
            if (!f){
                std::cout<<"Can't open "<<fname<<" for writing!\n";
            }
            for (int i=0; i<lat_size;i++){
                double v = cellLattice[i].vx();
                std::string endChar = (i%l==(l-1))? "\n":"\t";
                f<<v<<endChar;
            }
            f.close();
        }
        void latticeReset(){
            #pragma omp parallel for
            for (int j = 0; j<lat_size;j++){
                cellLattice[j].cellReset();
            }             
        }       
};
}
int main(int argc, char* argv[]){
    
    double u0=0.2;
    double alpha=0.02;
    double tau = 3.0*alpha+0.5; 

    if (argc!=6){
        std::cout<<"Hey man, 5 cml arguments\n";
        std::cout<<"<.exe><maxnThreads><length><width><timeStepsNumber><saveFlag>\n";
        return 0;
    }

    int maxnThreads = std::stoi(argv[1]);
    int l = std::stoi(argv[2]);
    int w = std::stoi(argv[3]);
    int time = std::stoi(argv[4]);
    int saveFlag = std::stoi(argv[5]);
    double *results = new double[maxnThreads];
    LBM::LatticeBoltzmann *mySim = new LBM::LatticeBoltzmann(l, w, tau); 

    for (int i = 1; i< maxnThreads+1;i++){
        results[i-1]=mySim->simulate(u0,time,i);
        if (i == saveFlag){
        // just one set of output for a random threadNum for sanity check
        std::string inf3 = "rho_02u0_"+std::to_string(l)+"x"+std::to_string(w)+"_"+std::to_string(time)+"sec_"+std::to_string(i)+"threads.txt";
        std::string inf4 = "xVel_02u0_"+std::to_string(l)+"x"+std::to_string(w)+"_"+std::to_string(time)+"sec_"+std::to_string(i)+"threads.txt";
        mySim->exportRho(inf3);
        mySim->exportxVel(inf4);
        }
        mySim->latticeReset();
        std::cout<<"thredNum = "<<i<<". Took "+std::to_string(results[i-1])+" seconds\n";
    }
    // save results to .txt file
    std::string fname = "LBM_"+std::to_string(l)+"x"+std::to_string(w)+"_"+std::to_string(time)+"steps_vs_nThreads.txt";
    std:: ofstream file(fname);
    if (!file.is_open()){
        printf("Error! Can't open file!\n");
        return 0;
    }
    else{
        file <<"nThreads\tTime\n";
        for (int i=0;i<maxnThreads;i++){
            file<<(i+1)<<"\t"<<results[i]<<"\n";
        } 
    }
    file.close();
    printf("Wrote results to file successfuly!\n");
    std::cout<<"Done!\n";
    return 1;
}
