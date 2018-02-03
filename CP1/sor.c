#include "cp1headers.h"

void sor(int nx, int ny, int niter, float* aw, float* ae, float* an, float* as, float* ap, float* p, float* s, float* ph, float* ttemp, float omega){
	
	float term;
	int i,j,ij,ijw,ije,ijn,ijs,iter;
	
	for (iter = 0; iter < niter; iter++){
		for ( i = 1; i < nx + 1; i++){
			for ( j = 1; j < ny + 1; j++){
			  ij  = i +j*(nx+2);
			  ije = (i+1)+ j*(nx+2);
			  ijw = (i-1)+ j*(nx+2);
			  ijs = i + (j-1)*(nx+2);
			  ijn = i + (j+1)*(nx+2);
			  term = aw[ij]*p[ijw] + ae[ij]*p[ije] + as[ij]*p[ijs] + an[ij]*p[ijn] + s[ij];
			  ph[ij] = omega * term /ap[ij] + (1.0 - omega) *p[ij];
			}
		}
		
		//Applying the Boundary condition for the pressure equation
		for ( j = 1; j < ny + 1; j++){
		  i = 0;
		  ij  = i +j*(nx+2); ije = (i+1)+ j*(nx+2);
		  ph[ij] = ph[ije];

		  i = nx+1;
		  ij  = i +j*(nx+2); ijw = (i-1)+ j*(nx+2);
		  ph[ij] = ph[ijw];
		}
		
		for ( i = 1; i < nx + 1; i++){
		  j = 0;
		  ij  = i +j*(nx+2); ijn = i + (j+1)*(nx+2);
		  ph[ij] = ph[ijn];

		  j = ny+1;
		  ij  = i +j*(nx+2);ijs = i + (j-1)*(nx+2);
		  ph[ij] = ph[ijs];
		}
		
	//Swap p with ph
    ttemp = ph; ph = p; p = ttemp;
    }
}