#pragma once
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
        static std::array<double, 9>f_i;
    public:
        Cell(){
            f_i = {0,0,0,0,0,0,0,0,0};
        }
        ~Cell();
        inline void cellInit(){
            // have some initial flow to the right
            f_i = {0,1,0,0,0,0,0,0,0};
        }
        // return value of macro density for a cell
        inline double density(){
            double ret = 0;
            // small loops like these are not worth parallelising
            for (int i=0; i<f_i.size();i++){
                ret+= f_i[i];
            }
            return ret;
        }
        // return value of macro flow velocity out of a cell
        inline std::pair<double, double> velocity(){
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
        inline std::array<double, 9> equilibrium(){
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
        // return all terms of the BGK collision operator of the cell
        inline std::array<double,9> BGK(double t){
            std::array<double,9> BGK;
            const std::array<double, 9> fiEQ=equilibrium();
            for (int i=0;i<f_i.size();i++){
                BGK[i] = -(f_i[i]-fiEQ[i])/t;
            }
            return std::move(BGK);
        }
};
}