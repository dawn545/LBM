#include "LBM.h"
#include <cstring>
#include <algorithm>
#include <cmath>

LBM::LBM(int N, float tau) : N(N), tau(tau), u0(0.06f) {
    int total_cells = N * N;
    // 分配内存
    f = new double[total_cells * 9];
    f_new = new double[total_cells * 9];
    rho = new float[total_cells];
    u = new float[total_cells];
    v = new float[total_cells];
    obstacle = new bool[total_cells];
    display_scalar = new float[total_cells];

    // D2Q9 模型标准参数定义
    // 方向约定：0为静止，1-4为上下左右，5-8为四个对角线
    double w_init[] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    int dir_x[] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    int dir_y[] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
    // 定义反方向数组：例如向右(1)的相反是向左(3)，主要用于固体表面的“反弹”处理
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
    // 释放内存，防止内存泄漏
    delete[] f; delete[] f_new;
    delete[] rho; delete[] u; delete[] v;
    delete[] obstacle; delete[] display_scalar;
}

void LBM::init() {
    // 遍历整个网格域进行初始化
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            int g_idx = grid_idx(x, y);
            rho[g_idx] = 1.0f;     // 初始密度设为 1.0
            u[g_idx] = u0;         // 初始具有向右的均匀流速 u0
            v[g_idx] = 0.0f;       // Y方向流速为 0
            obstacle[g_idx] = false;

            // 根据初始宏观速度和密度，计算每个方向的初始平衡态分布函数 (f_eq)
            for (int i = 0; i < 9; ++i) {
                // cu 代表速度向量点乘：3.0 * (e_i · u)
                double cu = 3.0 * (cx[i] * u[g_idx] + cy[i] * v[g_idx]);
                // 标准 LBM 的平衡态公式 (Taylor 展开后的 Maxwell-Boltzmann 分布)
                // feq = w_i * rho * [1 + 3(e_i·u) + 4.5(e_i·u)^2 - 1.5|u|^2]
                f[idx(x, y, i)] = w[i] * 1.0 * (1.0 + cu + 0.5 * cu * cu - 1.5 * (u[g_idx] * u[g_idx]));
                f_new[idx(x, y, i)] = f[idx(x, y, i)];
            }
        }
    }
}

void LBM::setObstacle() {
    // 在流场靠左的位置设置一个圆形障碍物，用于激发卡门涡街
    int cx_obs = N / 4;        // 圆心 X 坐标
    int cy_obs = N / 2;        // 圆心 Y 坐标
    int radius = N / 12;       // 圆柱半径
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            // 用圆的方程判定网格是否在障碍物内部
            if ((x - cx_obs) * (x - cx_obs) + (y - cy_obs) * (y - cy_obs) <= radius * radius) {
                obstacle[grid_idx(x, y)] = true;
            }
        }
    }
}

void LBM::step() {
    // 一个标准 LBM 时间步的三个流水线动作
    collideAndStream();         // 1. 碰撞与流迁
    applyBoundaries();          // 2. 处理计算域边界
    updateMacroscopicAndVis();  // 3. 更新宏观物理量，为渲染做准备
}

void LBM::collideAndStream() {
    // 遍历所有网格
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            int g_idx = grid_idx(x, y);
            // 固体内部不参与流体计算
            if (obstacle[g_idx]) continue;

            float r_c = rho[g_idx];
            float u_c = u[g_idx];
            float v_c = v[g_idx];

            for (int i = 0; i < 9; ++i) {
                // --- 步骤 1：计算碰撞 (Collide) ---
                // 基于局部宏观量，计算当前格子的平衡态 feq
                double cu = 3.0 * (cx[i] * u_c + cy[i] * v_c);
                double feq = w[i] * r_c * (1.0 + cu + 0.5 * cu * cu - 1.5 * (u_c * u_c + v_c * v_c));
                
                // BGK 近似：粒子分布函数 f 以 tau 为松弛时间，向平衡态 feq 靠拢
                // f_post 是碰撞后、流迁前的中间态
                double f_post = f[idx(x, y, i)] - (f[idx(x, y, i)] - feq) / tau;

                // --- 步骤 2：流迁 (Stream) ---
                // 计算粒子按照方向 i 移动一步后，将到达的邻居坐标 (nx, ny)
                int nx = x + cx[i];
                int ny = y + cy[i];
                
                // 处理 Y 方向的周期性边界条件 (上下连通，像吃豆人地图)
                if (ny < 0) ny = N - 1; 
                if (ny >= N) ny = 0;

                // 处理流体内部及障碍物边界
                if (nx >= 0 && nx < N) {
                    if (obstacle[grid_idx(nx, ny)]) {
                        // 【标准反弹边界 (Bounce-back)】
                        // 如果邻居是固体障碍物，粒子撞上去后原路弹回！
                        // 所以原本想流向 i 方向的粒子，最终流向了当前格子的 opposite[i] 方向
                        f_new[idx(x, y, opposite[i])] = f_post;
                    } else {
                        // 正常流迁：粒子成功抵达邻居格子 (nx, ny) 的 i 方向槽位中
                        f_new[idx(nx, ny, i)] = f_post;
                    }
                }
            }
        }
    }
    // 将更新好的 f_new 复制回 f，为下一帧做准备
    std::memcpy(f, f_new, sizeof(double) * N * N * 9);
}

void LBM::applyBoundaries() {
    // 处理左边界 (x=0)：强制施加一个稳定向右的速度 u0 (Zou-He 速度边界条件的简化版)
    for (int y = 0; y < N; ++y) {
        for (int i = 0; i < 9; ++i) {
            double cu = 3.0 * (cx[i] * u0);
            // 强行把左边界的分布函数重置为速度为 u0 的平衡态
            f[idx(0, y, i)] = w[i] * 1.0 * (1.0 + cu + 0.5 * cu * cu - 1.5 * (u0 * u0));
            
            // 处理右边界 (x=N-1)：简单的零梯度（流出）边界条件
            // 假设右边界的分布函数等于它左边一列的分布函数
            f[idx(N - 1, y, i)] = f[idx(N - 2, y, i)];
        }
    }
}

void LBM::updateMacroscopicAndVis() {
    // 遍历网格，从微观分布函数计算宏观流体力学量
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            int g_idx = grid_idx(x, y);
            
            // 障碍物内部宏观量归零
            if (obstacle[g_idx]) {
                rho[g_idx] = 1.0f; u[g_idx] = 0; v[g_idx] = 0;
                display_scalar[g_idx] = 0.0f; // 渲染为黑色
                continue;
            }
            
            double r = 0, ux = 0, uy = 0;
            // 对 9 个方向的分布函数进行统计积分（实际上就是离散求和）
            for (int i = 0; i < 9; ++i) {
                double val = f[idx(x, y, i)];
                r += val;                    // 零阶矩：总密度
                ux += cx[i] * val;           // 一阶矩的 X 分量：X方向动量
                uy += cy[i] * val;           // 一阶矩的 Y 分量：Y方向动量
            }
            
            rho[g_idx] = (float)r;
            // 速度 = 动量 / 密度
            u[g_idx] = (float)(ux / r);
            v[g_idx] = (float)(uy / r);
            
            // 计算当前格子的绝对速度大小
            float speed = std::sqrt(u[g_idx]*u[g_idx] + v[g_idx]*v[g_idx]);
            // 将速度归一化到 0.0 ~ 1.0 之间，存入 display_scalar 供 OpenGL 渲染颜色
            display_scalar[g_idx] = std::fmin(1.0f, speed / u0);
        }
    }
}
