#pragma once
#include <array>
class Cell{
    private:
        const std::array<double, 9> weights {{4./9., 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36}};
        // like the diagram in the lecture notes: static, E, W, N, S, NE, SW, NW, SE
        static std::array<std::pair<int,int>, 9> e_i {{0,0},{1,0},{-1,0},{0,1},{0,-1},{1,1},{-1,-1},{-1,1},{1,-1}}
    public:
        static::array<double, 9>f_i;
        Cell();
        ~Cell();
        double density();
        double std::pair<double, double> velocity();
        std::array<double, 9> equilibrium(double accelXtau);
        void collision(double accelXtau, double tau);
}