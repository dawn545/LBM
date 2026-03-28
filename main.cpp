#include <GL/freeglut.h>
#include <iostream>
#include "LBM.h"

const int N = 128;
const int STEPS_PER_FRAME = 10;

// 修复语法：类型 变量名(构造参数);
LBM fluid(N, 0.53f); 

void init() {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    float h = 1.0f / N;
    float* scalar = fluid.getDensity();

    glBegin(GL_QUADS);
    for (int y = 0; y < N - 1; y++) {
        for (int x = 0; x < N - 1; x++) {
            float d00 = scalar[y * N + x];
            float d10 = scalar[y * N + x + 1];
            float d11 = scalar[(y + 1) * N + x + 1];
            float d01 = scalar[(y + 1) * N + x];

            float px = x * h, py = y * h;
            glColor3f(d00 * 0.8f, d00 * 0.9f, d00 + 0.2f); glVertex2f(px, py);
            glColor3f(d10 * 0.8f, d10 * 0.9f, d10 + 0.2f); glVertex2f(px + h, py);
            glColor3f(d11 * 0.8f, d11 * 0.9f, d11 + 0.2f); glVertex2f(px + h, py + h);
            glColor3f(d01 * 0.8f, d01 * 0.9f, d01 + 0.2f); glVertex2f(px, py + h);
        }
    }
    glEnd();
    glutSwapBuffers();
}

void idle() {
    for (int i = 0; i < STEPS_PER_FRAME; i++) fluid.step();
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    if (key == 27 || key == 'q') exit(0);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutCreateWindow("LBM Flow");
    init();
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}