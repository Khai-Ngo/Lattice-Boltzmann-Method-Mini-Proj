#include <vector>
#include <array>
#include <iostream>
#include <utility>
#include <fstream>
#include <cmath>
#include <string>

namespace LBM{
class Cell{
    private:
        const std::array<double, 9> weights {{4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36}};
        // like the diagram in the lecture notes: static, E, W, N, S, NE, SW, NW, SE
        const std::array<int, 9>ex {{0,1,-1,0,0,1,-1,-1,1}};
        const std::array<int, 9>ey {{0,0,0,1,-1,1,-1,1,-1}};
        std::array<double, 9>f_i;
        std::array<double, 9>f_dup; // for updating. Avoid overwriting data before they get propagated
    public:
        Cell(){
             // concurrent f_i and its duplicate. Default to same values as the weights (recommended by textbooks)
            this->f_i = {4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};
            this->f_dup = {4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36};
        }
        // if you want a different initial flow, specify it yourself
        Cell(const std::array<double, 9>& f_i){
            this->f_i = f_i;
            this->f_dup = f_i;
        }
        Cell(double rho){
            // init with zero velocity. Basically the weights array times a multiplier
            this-> f_i={rho*4./9., rho*1./9, rho*1./9, rho*1./9, rho*1./9, rho*1./36, rho*1./36, rho*1./36, rho*1./36};
            this-> f_dup={rho*4./9., rho*1./9, rho*1./9, rho*1./9, rho*1./9, rho*1./36, rho*1./36, rho*1./36, rho*1./36};
        }
        // copy constructor
        Cell (const Cell& inputCell){
            this->f_i = inputCell.f_i;
            this->f_dup = inputCell.f_dup;
        }
        void update_fi(int dir, double input){
            this->f_dup[dir] = input;
        }
        double get_fi(int dir){
            return this->f_i[dir];
        }
        void sync(){
            this->f_i = this->f_dup;
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
        // return macro speed at the cell
        double speed(){
            double vx=0;
            double vy=0;
            double rho=0;
            for (int i=0;i<f_i.size();i++){
                vx+=ex[i]*f_i[i];
                vy+=ey[i]*f_i[i];
                rho+=f_i[i]; 
            }
            vx /= rho;
            vy /= rho;
            return std::sqrt(vx*vx+vy*vy);
        }
        // Avoid looping through the array many times but repeat code
        void collision(double t){
            // don't call speed() or velocity(), because if you do you need call density() for rho, and that means looping through f_i twice
            double vx=0;
            double vy=0;
            double rho=0;
            for (int i=0;i<f_i.size();i++){
                vx+=ex[i]*f_i[i];
                vy+=ey[i]*f_i[i];
                rho+=f_i[i]; 
            }
            vx /= rho;
            vy /= rho;
            double dotProd, uProd, fiEQ_i;
            for (int i = 0; i<f_i.size();i++){  
                dotProd = ex[i]*vx+ey[i]*vy;
                uProd = vx*vx+vy*vy;
                fiEQ_i= weights[i]*rho*(1+3*dotProd+4.5*dotProd*dotProd-1.5*uProd);
                f_i[i] -= (f_i[i]-fiEQ_i)/t;
            }
        }
        // for debugging purposes
        std::array<double, 9> distFuctions(){
            return f_i;
        }
        // return value of macro flow velocity out of a cell
        std::pair<double, double> velocity(){
            double pvx=0;
            double pvy=0;
            double rho=0;
            for (int i=0;i<f_i.size();i++){
                pvx+=ex[i]*f_i[i];
                pvy+=ey[i]*f_i[i];
                rho+=f_i[i]; // dont call density() to avoid looping it twice
            }
            return std::make_pair (pvx/rho, pvy/rho); 
        }
        std::array<double, 9> equilibrium(){
            std::array<double, 9> fiEQ;
            std::pair v = velocity();
            double vx = std::get<0>(v);
            double vy = std::get<1>(v);
            double rho = density(); // oh well, for cleaniness's sake
            for (int i = 0; i<f_i.size();i++){  
                double dotProd = ex[i]*vx+ey[i]*vy;
                double uProd = vx*vx+vy*vy;
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
        int nodeType(int i){
            // find out what type of node it is based on node index
            if (isTop(i)){
                        if (isLeft(i)){
                            return 1; // top left node
                        }
                        else if (isRight(i))
                        {
                            return 2; // top right node
                        }
                        else return 3; // top row node
                        
                    }
                    else if (isBottom(i))
                    {
                        if (isLeft(i)){
                            return 4; // bot left node
                        }
                        else if (isRight(i))
                        {
                            return 5; // bot right node
                        }
                        else return 6; // bot row node
                    }
                    else if (isLeft(i)){
                        return 7; // left col node
                    }
                    else if (isRight(i))
                    {
                        return 8; // right col node 
                    }
                    else return 9; // middle node
                }
        // 9 (kinda) types of nodes with different propagation rules (due to them being affected by BC in different ways)
        void propTopLeft(int i){
            this->cellLattice[i+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
            this->cellLattice[i-1+l].update_fi(2, cellLattice[i].get_fi(2)); // e2 W affected by periodic BC
            this->cellLattice[i].update_fi(4, cellLattice[i].get_fi(3)); // e3 N reflect off top wall
            this->cellLattice[i+l].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
            this->cellLattice[i+2].update_fi(8, cellLattice[i].get_fi(5)); // e5 NE specular reflect off top wall
            this->cellLattice[i+2*l-1].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW affected by periodic BC
            this->cellLattice[i+l-2].update_fi(6, cellLattice[i].get_fi(7)); // e7 NW affected by periodic BC and specular reflect off top wall
            this->cellLattice[i+l+1].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE 
        }
        void propTopRight(int i){
            this->cellLattice[i-l+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E affected by periodic BC
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
            this->cellLattice[i].update_fi(4, cellLattice[i].get_fi(3)); // e3 N reflect off top wall
            this->cellLattice[i+l].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
            this->cellLattice[i-l+2].update_fi(8, cellLattice[i].get_fi(5)); // e5 NE affected by periodic BC and specular reflect off top wall
            this->cellLattice[i+l-1].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW
            this->cellLattice[i-2].update_fi(6, cellLattice[i].get_fi(7)); // e7 NW specular reflect off top wall
            this->cellLattice[i+1].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE PBC
        }
        void propBotLeft(int i){
            this->cellLattice[i+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W affected by periodic BC 
            this->cellLattice[i-l].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
            this->cellLattice[i].update_fi(3, cellLattice[i].get_fi(4)); // e4 S reflect off bot wall
            this->cellLattice[i-l+1].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE
            this->cellLattice[i+l-2].update_fi(7, cellLattice[i].get_fi(6)); // e6 SW specular refl off bot wall + PBC
            this->cellLattice[i-1].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW PBC
            this->cellLattice[i+2].update_fi(5, cellLattice[i].get_fi(8)); // e8 SE specular refl off bot wall
        }
        void propBotRight(int i){
            this->cellLattice[i-l+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E PBC
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
            this->cellLattice[i-l].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
            this->cellLattice[i].update_fi(3, cellLattice[i].get_fi(4)); // e4 S reflect rev bot wall
            this->cellLattice[i-2*l+1].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE PBC
            this->cellLattice[i-2].update_fi(7, cellLattice[i].get_fi(6)); // e6 SW spec refl bot wall
            this->cellLattice[i-l-1].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW
            this->cellLattice[i-l+2].update_fi(5, cellLattice[i].get_fi(8)); // e8 SE PBC + spec refl bot wall 
        }
        void propLeftCol(int i){
            this->cellLattice[i+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
            this->cellLattice[i-1+l].update_fi(2, cellLattice[i].get_fi(2)); // e2 W PBC
            this->cellLattice[i-l].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
            this->cellLattice[i+l].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
            this->cellLattice[i-l+1].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE
            this->cellLattice[i+2*l-1].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW PBC
            this->cellLattice[i-1].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW PBC
            this->cellLattice[i+l+1].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE 
        }
        void propRightCol(int i){
            this->cellLattice[i-l+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E PBC
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
            this->cellLattice[i-l].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
            this->cellLattice[i+l].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
            this->cellLattice[i-2*l+1].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE PBC
            this->cellLattice[i+l-1].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW
            this->cellLattice[i-l-1].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW
            this->cellLattice[i+1].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE PBC
        }
        void propTopRow(int i){
            this->cellLattice[i+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
            this->cellLattice[i].update_fi(4, cellLattice[i].get_fi(3)); // e3 N refl rev top wall
            this->cellLattice[i+l].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
            int j = ((i%l + 2) >= l)? i+2-l:i+2; // pesky 2-unit shift by spec refl...
            this->cellLattice[j].update_fi(8, cellLattice[i].get_fi(5)); // e5 NE spec refl top wall
            this->cellLattice[i+l-1].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW 
            int k = ((i%l - 2) < 0)? i-2+l:i-2; 
            this->cellLattice[k].update_fi(6, cellLattice[i].get_fi(7)); // e7 NW spec refl top wall
            this->cellLattice[i+l+1].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE
        }
        void propBotRow(int i){
            this->cellLattice[i+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
            this->cellLattice[i-l].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
            this->cellLattice[i].update_fi(3, cellLattice[i].get_fi(4)); // e4 S rfl rev bot wall
            this->cellLattice[i-l+1].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE
            int j = ((i%l-2)<0)? i-2+l:i-2;
            this->cellLattice[j].update_fi(7, cellLattice[i].get_fi(6)); // e6 SW spec refl bot wall
            this->cellLattice[i-l-1].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW
            int k = ((i%l + 2) >= l)? i+2-l:i+2;
            this->cellLattice[k].update_fi(5, cellLattice[i].get_fi(8)); // e8 SE spec refl bot wall
        }
        void propMid(int i){
            this->cellLattice[i+1].update_fi(1, cellLattice[i].get_fi(1)); // e1 E
            this->cellLattice[i-1].update_fi(2, cellLattice[i].get_fi(2)); // e2 W 
            this->cellLattice[i-l].update_fi(3, cellLattice[i].get_fi(3)); // e3 N
            this->cellLattice[i+l].update_fi(4, cellLattice[i].get_fi(4)); // e4 S
            this->cellLattice[i-l+1].update_fi(5, cellLattice[i].get_fi(5)); // e5 NE
            this->cellLattice[i+l-1].update_fi(6, cellLattice[i].get_fi(6)); // e6 SW
            this->cellLattice[i-l-1].update_fi(7, cellLattice[i].get_fi(7)); // e7 NW
            this->cellLattice[i+l+1].update_fi(8, cellLattice[i].get_fi(8)); // e8 SE
        }
    public:
        LatticeBoltzmann(){};
        LatticeBoltzmann(int l, int w, double tau){
            this->l = l;
            this->w = w;
            this->tau = tau;
            this->lat_size = l*w;
            // initialize all lattice cells to default distributions. Parallelisable
            Cell *initCell = new Cell(); // does not look neat. Leave it for now
            for (int i=0; i<lat_size; i++){        
                this->cellLattice.push_back(*initCell);
            }
            delete initCell;      
        }
        void initctrPressurePulse(int l, int w, double tau){
            this->l = l;
            this->w = w;
            this->tau = tau;
            this->lat_size = l*w;
            int ctr_index = xy_to_index(l/2, w/2, l);
            for (int i=0; i<lat_size;i++){
                double rho = static_cast<float>((i/l+1)*(l-i/l)*(i%l+1)*(w-i%l))/static_cast<float>(lat_size*lat_size/4);
                Cell *temp = new Cell(rho);
                this->cellLattice.push_back(*temp);
                delete temp;
            }
        }
        void simulate(int time){
            for (int t=0; t<time;t++){
                for (int i=0;i<lat_size;i++){
                    // collision step
                    this->cellLattice[i].collision(tau);
                    // advection step
                    switch (nodeType(i))
                    {
                    case 1:
                        propTopLeft(i);
                        break;
                    case 2:
                        propTopRight(i);
                        break;
                    case 3:
                        propTopRow(i);
                        break;
                    case 4:
                        propBotLeft(i);
                        break;
                    case 5:
                        propBotRight(i);
                        break;
                    case 6:
                        propBotRow(i);
                        break;
                    case 7:
                        propLeftCol(i);
                        break;
                    case 8:
                        propRightCol(i);
                        break;
                    case 9:
                        propMid(i);
                        break;
                    default:
                        std::cout<<"You should never see this line printed\n";
                        break;
                    } 
                }
                // sync current flux with upated flux
                for (int j = 0; j<lat_size;j++)
                    cellLattice[j].sync();
            }
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
        void exportSpd(std::string& fname){
            std::ofstream f (fname);
            if (!f){
                std::cout<<"Can't open "<<fname<<" for writing!\n";
            }
            for (int i=0; i<lat_size;i++){
                double v = cellLattice[i].speed();
                std::string endChar = (i%l==(l-1))? "\n":"\t";
                f<<v<<endChar;
            }
            f.close();
        }
        
};

/*
std::pair<int, int> index_to_xy(int index, int l){
    return std::make_pair(index/l, index%l);
}
*/
}
int main(int argc, char* argv[]){
    int l = 10;
    int w = 10;
    double tau = 0.5;  
    int time = 1000;
    
    LBM::LatticeBoltzmann *mySim = new LBM::LatticeBoltzmann();
    mySim->initctrPressurePulse(l, w, tau);  
    
    std::string inf1("init_rho_ctrpuls.txt");
    std::string inf2("inti_spd_ctrpuls.txt");
    mySim->exportRho(inf1);
    mySim->exportSpd(inf2);
    
    mySim->simulate(time);
    std::string inf3("fin_rho_ctrpuls_10x10_050tau_1ks.txt");
    std::string inf4("fin_spd_ctrpuls_10x10_050tau_1ks.txt");
    mySim->exportRho(inf3);
    mySim->exportSpd(inf4);
    std::cout<<"Done!\n";
    return 0;
}
