#include "cp1headers.h"

void ppesource(int nx, int ny, float* uh, float* vh, float* s, float dx, float dy, float dt){
	
	int i,j,ij,ijw,ije,ijn,ijs;
	
	for ( i = 1; i < nx + 1; i++){
		for ( j = 1; j < ny + 1; j++){
			ij  = i +j*(nx+2);
			ije = (i+1)+ j*(nx+2);
			ijw = (i-1)+ j*(nx+2);
			ijs = i + (j-1)*(nx+2);
			ijn = i + (j+1)*(nx+2);
			s[ij] = (uh[ijw] - uh[ije])/(2.0*dx) + (vh[ijs] - vh[ijn])/(2.0*dy);
			s[ij] = s[ij]/dt;
		}
    }
	
}