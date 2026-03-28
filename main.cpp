#include <GL/freeglut.h>
#include <iostream>
#include "LBM.h"

const int N = 128;               // 模拟的网格分辨率 128x128
// 极其重要：时间加速器。LBM 时间步非常小，如果每算一步就渲染一帧，画面会像幻灯片一样慢。
// 这里设置为每计算 10 步流体演化，屏幕才更新 1 帧。
const int STEPS_PER_FRAME = 10;  

// 全局实例化 LBM 模拟器。tau = 0.53 意味着极低的粘度 (运动学粘度 nu = (tau-0.5)/3 = 0.01)
// 这很容易产生卡门涡街。
LBM fluid(N, 0.53f); 

void init() {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // 背景涂黑
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);       // 设置正交投影，将坐标系映射到 0.0~1.0 之间
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    float h = 1.0f / N; // 每个网格在屏幕上的长宽比例因子
    
    // 拿到我们在 updateMacroscopicAndVis() 里算好的 0~1 的速度标量场
    float* scalar = fluid.getDensity();

    // 开启四边形绘制模式
    glBegin(GL_QUADS);
    // 遍历整个网格（留出边缘1格防止越界）
    for (int y = 0; y < N - 1; y++) {
        for (int x = 0; x < N - 1; x++) {
            // 获取当前网格及周围三个网格的标量值 (用于平滑着色)
            float d00 = scalar[y * N + x];
            float d10 = scalar[y * N + x + 1];
            float d11 = scalar[(y + 1) * N + x + 1];
            float d01 = scalar[(y + 1) * N + x];

            // 映射到屏幕坐标系 (0.0 ~ 1.0)
            float px = x * h, py = y * h;
            
            // 基于标量值设置颜色，这里用了一个简单的调色板映射：
            // RGB = (val*0.8, val*0.9, val+0.2) 
            // 速度(val)越趋近 0，颜色越深蓝；速度越趋近 1，颜色越偏向亮青/白色。
            glColor3f(d00 * 0.8f, d00 * 0.9f, d00 + 0.2f); glVertex2f(px, py);
            glColor3f(d10 * 0.8f, d10 * 0.9f, d10 + 0.2f); glVertex2f(px + h, py);
            glColor3f(d11 * 0.8f, d11 * 0.9f, d11 + 0.2f); glVertex2f(px + h, py + h);
            glColor3f(d01 * 0.8f, d01 * 0.9f, d01 + 0.2f); glVertex2f(px, py + h);
        }
    }
    glEnd();
    
    // 交换前后缓冲区，将绘制好的画面显示到屏幕上
    glutSwapBuffers();
}

void idle() {
    // OpenGL 的空闲回调函数（只要 CPU 没事干就会不停执行）
    // 连续计算 N 个 LBM 时间步
    for (int i = 0; i < STEPS_PER_FRAME; i++) {
        fluid.step();
    }
    // 通知 OpenGL 重新调用 display() 绘制画面
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    // 按 Esc 或 Q 键退出程序
    if (key == 27 || key == 'q') exit(0);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); // 启用双缓冲和 RGB 模式
    glutInitWindowSize(800, 800);                // 渲染窗口分辨率 800x800
    glutCreateWindow("LBM Flow");                // 窗口标题
    
    init();
    
    // 注册回调函数
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    
    // 启动 GLUT 事件处理主循环
    glutMainLoop();
    return 0;
}
