#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void initialize(int nx, int ny, float* u, float* v, float* p);
void ppecoeffs(int nx, int ny, float* aw, float* ae, float* an, float* as, float* ap, float dx, float dy);
void momentum(int nx, int ny, float* u, float* v, float* uh, float* vh, float* s, float amu, float dx, float dy, float dt);
void bc(int nx, int ny, float* u, float* v, float utop);
void sor(int nx, int ny, int niter, float* aw, float* ae, float* an, float* as, float* ap, float* p, float* s, float* ph, float* ttemp, float omega);
void momentumcorr(int nx, int ny, float* u, float* v, float* p, float* uh, float* vh, float dx, float dy, float dt);
void ppesource(int nx, int ny, float* uh, float* vh, float* s, float dx, float dy, float dt);
void writeflowfield(int nx, int ny, float dx, float dy, float* u, float* v, float* p);
void writeparticlepos(float* xp, float* yp, float* up, float* vp, int step, int np);
void bilinterp(float *u, float *v, float *xp, float *yp, float* up, float *vp, int np, int nx, int ny, float dx, float dy);