#include "cp1headers.h"

void initialize(int nx, int ny, int np, float* u, float* v, float* p, float* up, float* vp){
	int i,j,ij;
	
	for ( i = 0; i < nx + 2; i++){
		for ( j = 0; j < ny + 2; j++){
		  ij = i +j*(nx+2);
		  u[ij] = 0.0;
		  v[ij] = 0.0;
		  p[ij] = 0.0;
		}
	}
	
	for (i = 0; i<np; i++){
		up[i] = 0.0;
		vp[i] = 0.0;
	}
	
}