#include "cp1headers.h"

void momentumcorr(int nx, int ny, float* u, float* v, float* p, float* uh, float* vh, float dx, float dy, float dt){
	
	int i,j,ij,ijw,ije,ijn,ijs;
	
	//Correct u and v for p	
	for ( i = 1; i < nx + 1; i++){
		for ( j = 1; j < ny + 1; j++){
			ij  = i +j*(nx+2);
			ije = (i+1)+ j*(nx+2);
			ijw = (i-1)+ j*(nx+2);
			ijs = i + (j-1)*(nx+2);
			ijn = i + (j+1)*(nx+2);
			u[ij] = uh[ij] + dt * (p[ijw] - p[ije])/(2.0*dx);
			v[ij] = vh[ij] + dt * (p[ijs] - p[ijn])/(2.0*dy);
		}
    }
	
}