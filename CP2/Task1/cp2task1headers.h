#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void initialize(int nx, int ny, float dx, float dy, float xcenter0, float ycenter0, float radius, float velocity, float* u, float* v, float* phi);
void bc(int nx, int ny, float* u, float* v, float* phi);
void writeflowfield(int nx, int ny, float dx, float dy, float* u, float* v);
void writephi(int nx, int ny, float dx, float dy, float* phi, int step);
void swap(int nx, int ny, float* phi, float* phi_new);
void advect(int nx, int ny, float dx, float dy, float dt, float* u, float* v, float* phi, float* phi_new);
float sign(float phi);