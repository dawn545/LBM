#include "LBM.h"
#include <cstring>
#include <algorithm>
#include <cmath>

LBM::LBM(int N, float tau) : N(N), tau(tau), u0(0.06f) {
    int total_cells = N * N;
    f = new double[total_cells * 9];
    f_new = new double[total_cells * 9];
    rho = new float[total_cells];
    u = new float[total_cells];
    v = new float[total_cells];
    obstacle = new bool[total_cells];
    display_scalar = new float[total_cells]; // 修复拼写

    double w_init[] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    int dir_x[] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    int dir_y[] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
    int opp[]   = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

    for (int i = 0; i < 9; ++i) {
        w[i] = w_init[i];
        cx[i] = dir_x[i];
        cy[i] = dir_y[i];
        opposite[i] = opp[i];
    }
    init();
    setObstacle();
}

LBM::~LBM() {
    delete[] f; delete[] f_new;
    delete[] rho; delete[] u; delete[] v;
    delete[] obstacle; delete[] display_scalar;
}

void LBM::init() {
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            int g_idx = grid_idx(x, y);
            rho[g_idx] = 1.0f;
            u[g_idx] = u0;
            v[g_idx] = 0.0f;
            obstacle[g_idx] = false;
            for (int i = 0; i < 9; ++i) {
                double cu = 3.0 * (cx[i] * u[g_idx] + cy[i] * v[g_idx]);
                f[idx(x, y, i)] = w[i] * 1.0 * (1.0 + cu + 0.5 * cu * cu - 1.5 * (u[g_idx] * u[g_idx]));
                f_new[idx(x, y, i)] = f[idx(x, y, i)];
            }
        }
    }
}

void LBM::setObstacle() {
    int cx_obs = N / 4;
    int cy_obs = N / 2;
    int radius = N / 12;
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            if ((x - cx_obs) * (x - cx_obs) + (y - cy_obs) * (y - cy_obs) <= radius * radius) {
                obstacle[grid_idx(x, y)] = true;
            }
        }
    }
}

void LBM::step() {
    collideAndStream();
    applyBoundaries();
    updateMacroscopicAndVis();
}

void LBM::collideAndStream() {
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            int g_idx = grid_idx(x, y);
            if (obstacle[g_idx]) continue;

            float r_c = rho[g_idx];
            float u_c = u[g_idx];
            float v_c = v[g_idx];

            for (int i = 0; i < 9; ++i) {
                double cu = 3.0 * (cx[i] * u_c + cy[i] * v_c);
                double feq = w[i] * r_c * (1.0 + cu + 0.5 * cu * cu - 1.5 * (u_c * u_c + v_c * v_c));
                double f_post = f[idx(x, y, i)] - (f[idx(x, y, i)] - feq) / tau;

                int nx = x + cx[i];
                int ny = y + cy[i];
                if (ny < 0) ny = N - 1; if (ny >= N) ny = 0;

                if (nx >= 0 && nx < N) {
                    if (obstacle[grid_idx(nx, ny)]) {
                        f_new[idx(x, y, opposite[i])] = f_post;
                    } else {
                        f_new[idx(nx, ny, i)] = f_post;
                    }
                }
            }
        }
    }
    std::memcpy(f, f_new, sizeof(double) * N * N * 9);
}

void LBM::applyBoundaries() {
    for (int y = 0; y < N; ++y) {
        for (int i = 0; i < 9; ++i) {
            double cu = 3.0 * (cx[i] * u0);
            f[idx(0, y, i)] = w[i] * 1.0 * (1.0 + cu + 0.5 * cu * cu - 1.5 * (u0 * u0));
            f[idx(N - 1, y, i)] = f[idx(N - 2, y, i)];
        }
    }
}

void LBM::updateMacroscopicAndVis() {
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            int g_idx = grid_idx(x, y);
            if (obstacle[g_idx]) {
                rho[g_idx] = 1.0f; u[g_idx] = 0; v[g_idx] = 0;
                display_scalar[g_idx] = 0.0f;
                continue;
            }
            double r = 0, ux = 0, uy = 0;
            for (int i = 0; i < 9; ++i) {
                double val = f[idx(x, y, i)];
                r += val; ux += cx[i] * val; uy += cy[i] * val;
            }
            rho[g_idx] = (float)r;
            u[g_idx] = (float)(ux / r);
            v[g_idx] = (float)(uy / r);
            float speed = std::sqrt(u[g_idx]*u[g_idx] + v[g_idx]*v[g_idx]);
            display_scalar[g_idx] = std::fmin(1.0f, speed / u0);
        }
    }
}