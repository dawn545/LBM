#pragma once
#include <vector>
#include <cmath>

class LBM {
public:
    // 构造函数：N为网格分辨率(NxN)，tau为松弛时间(决定流体粘度)
    LBM(int N, float tau);
    ~LBM();

    // 执行一个完整的时间步（包含碰撞、流迁、边界处理、宏观量更新）
    void step();

    // Getter 方法，供外部（如 OpenGL 渲染端）获取模拟数据
    float* getDensity() { return display_scalar; } // 获取用于渲染的标量场（通常映射为速度大小或涡度）
    float* getVx() { return u; }                   // 获取 X 方向宏观速度场
    float* getVy() { return v; }                   // 获取 Y 方向宏观速度场
    bool* getObstacle() { return obstacle; }       // 获取障碍物掩码场

private:
    int N;        // 网格的边长（总网格数为 N * N）
    float tau;    // 松弛时间 (Relaxation time)。tau 越接近 0.5，流体粘度越低（越容易湍流）；tau 越大，粘度越高。
    float u0;     // 入口初始流速（也是系统特征速度）

    // --- LBM 核心数据数组 ---
    // 采用 D2Q9 模型（二维空间，每个格子有 9 个离散速度方向）
    double* f;      // 当前时刻的粒子分布函数 (Distribution functions)，大小为 N*N*9
    double* f_new;  // 下一时刻的粒子分布函数。为了防止在流迁(Streaming)时数据被覆盖，需要双缓冲
    
    // --- 宏观物理量 ---
    float* rho;     // 宏观密度场 (Density)，由 9 个 f 累加得到
    float* u;       // 宏观 X 速度场 (X-Velocity)
    float* v;       // 宏观 Y 速度场 (Y-Velocity)
    bool* obstacle; // 障碍物标记场 (true 表示该格子是固体边界)
    float* display_scalar; // 专用于图形渲染的标量场（本代码中存储了归一化后的速度大小）

    // --- D2Q9 模型常量 ---
    double w[9];        // 9 个方向的权重常数 (Weights)
    int cx[9];          // 9 个方向的离散速度 X 分量 (Velocity vectors X)
    int cy[9];          // 9 个方向的离散速度 Y 分量 (Velocity vectors Y)
    int opposite[9];    // 记录每个方向的相反方向索引（用于反弹边界条件 Bounce-back）

    // --- 辅助内联函数：将 2D/3D 坐标展平为 1D 数组索引 ---
    // 获取特定位置 (x,y) 上方向为 d 的分布函数 f 的一维索引
    inline int idx(int x, int y, int d) const { return (y * N + x) * 9 + d; }
    // 获取特定位置 (x,y) 的宏观量（如 rho, u, v）的一维索引
    inline int grid_idx(int x, int y) const { return y * N + x; }

    // --- LBM 算法核心步骤 ---
    void init();                     // 初始化流场（密度、速度及初始平衡态分布）
    void setObstacle();              // 放置物理障碍物（如圆柱）
    void collideAndStream();         // 核心：LBGK碰撞松弛与粒子相邻流迁
    void applyBoundaries();          // 处理计算域外部边界（如左侧固定流入）
    void updateMacroscopicAndVis();  // 积分分布函数计算新的宏观密度和速度
};
