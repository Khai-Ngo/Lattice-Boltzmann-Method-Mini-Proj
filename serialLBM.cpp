#include<vector>
#include <array>
#include <iostream>
#include <utility>

namespace LBM{
class Cell{
    private:
        const std::array<double, 9> weights {{4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36}};
        // like the diagram in the lecture notes: static, E, W, N, S, NE, NW, SW, SE
        const std::array<int, 9>ex {{0,1,-1,0,0,1,-1,-1,1}};
        const std::array<int, 9>ey {{0,0,0,1,-1,1,1,-1,-1}};
        std::array<double, 9>f_i;
    public:
        Cell(){
            // recommded way of initialising density from a Bsc diss. Means it's relatively static?
            f_i = weights;
        }
        Cell(std::array<double, 9> input){
            f_i = input;
        }
        // return value of macro density for a cell
        double density(){
            double ret = 0;
            // small loops like these are not worth parallelising
            for (int i=0; i<f_i.size();i++){
                ret+= f_i[i];
            }
            return ret;
        }
        // return value of macro flow velocity out of a cell
        std::pair<double, double> velocity(){
            double pvx=0;
            double pvy=0;
            double rho=0;
            for (int i=0;i<f_i.size();i++){
                pvx+=ex[i]*f_i[i];
                pvy+=ey[i]*f_i[i];
                rho+=f_i[i];
            }
            return std::make_pair (pvx/rho, pvy/rho); 
        }
        // return all f_i equilibrium terms of the cell, to be used for evolution
        std::array<double, 9> equilibrium(){
            std::array<double, 9> fiEQ;
            std::pair v = velocity();
            double vx = std::get<0>(v);
            double vy = std::get<1>(v);
            double rho = density();
            for (int i = 0; i<f_i.size();i++){  
                double dotProd = ex[i]*vx+ey[i]*vy;
                fiEQ[i]= weights[i]*rho*(1+3*dotProd+4.5*dotProd-1.5*dotProd);
            }
            return std::move(fiEQ); //  to avoid memory run-away
        }
        void collision(double t){
            std::array<double,9> BGK;
            const std::array<double, 9> fiEQ=equilibrium();
            for (int i=0;i<f_i.size();i++){
                BGK[i] = -(f_i[i]-fiEQ[i])/t;
                f_i[i] -= BGK[i];
            }
        }
        // debugging purposes
        std::array<double, 9> distFuctions(){
            return f_i;
        }
};


void readInputParams(std::string inf);
void printResults(std::vector<Cell> lattice, std::string outf);

std::pair<int, int> index_to_xy(int index, int l, int w){
    int x, y;
    x = index/l;
    y = index%l;
    return std::make_pair(x, y);
}

int main(int argc, char* argv[]){
    int grid_l = 200;
    int grid_w = 200;
    int lat_size = grid_l*grid_w;
    int time = 1000;
    int tau = 0.5;
    // initialize all lattice cells to default distributions. Parallelisable
    std::vector<Cell> Lattice;
    Cell tempCell; // does not look good but let's keep it for now
    for (int i =0; i<lat_size; i++){        
        Lattice.push_back(tempCell);
    }
    
    // create the obstacle

    // simulation main loop
    for (int t = 0; t<time;t++){
        // bulk computation here parallelisable
        for (int i=0;i<lat_size;i++){
            // collision step
            Lattice[i].collision(tau);
            // reflection step

            // advection step
        }
    }
    return 0;
}
}