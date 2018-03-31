#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void initialize(int nx, int ny, float dx, float dy, float xcenter0, float ycenter0, float radius, float velocity, float* u, float* v, float* phi);
void bcvelocity(int nx, int ny, float* u, float* v);
void bcphi(int nx, int ny, float* phi);
void swap(int nx, int ny, float* phi, float* phi_new);
void writeflowfield(int nx, int ny, float dx, float dy, float* u, float* v);
void writephi(int nx, int ny, float dx, float dy, float* phi, int step);
void advectx(int nx, int ny, float dx, float dt, float* u, float* phi, float* phi_new, float epsilon);
void advecty(int nx, int ny, float dy, float dt, float* v, float* phi, float* phi_new, float epsilon);
void reinitialize(int nx, int ny, int niter, float dx, float dy, float dtau, float* phi, float* phi_new, float* phi0, float epsilon);
float sign(float phi, float epsilon);