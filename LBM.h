#pragma once
#include <vector>
#include <cmath>

class LBM {
public:
    LBM(int N, float tau);
    ~LBM();

    void step();

    float* getDensity() { return display_scalar; }
    float* getVx() { return u; }
    float* getVy() { return v; }
    bool* getObstacle() { return obstacle; }

private:
    int N;
    float tau;
    float u0; 

    double* f;      
    double* f_new;  
    float* rho;     
    float* u;       
    float* v;       
    bool* obstacle; 
    float* display_scalar; // 确保这里拼写正确

    double w[9];
    int cx[9];
    int cy[9];
    int opposite[9];

    inline int idx(int x, int y, int d) const { return (y * N + x) * 9 + d; }
    inline int grid_idx(int x, int y) const { return y * N + x; }

    void init();
    void setObstacle();
    void collideAndStream();
    void applyBoundaries();
    void updateMacroscopicAndVis();
}; 