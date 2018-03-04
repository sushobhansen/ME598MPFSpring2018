#include "cp1headers.h"

void bilinterp(float *u, float *v, float *xp, float *yp, float* up, float *vp, int np, int nx, int ny, float dx, float dy){
	
	int i,j,ij,ije,ijn,ijne,k;
	
	for (k=0; k<np; k++){
		
		//Locate bottom left grid point for each particle
		i = (int)floor(xp[k]/dx);
		j = (int)floor(yp[k]/dy);
		
		ij  = i + j*(nx+2);
		ije = (i+1) + j*(nx+2);
		ijn = i + (j+1)*(nx+2);
		ijne = (i+1) + (j+1)*(nx+2);
		
		//Bilinear interpolation - see https://en.wikipedia.org/wiki/Bilinear_interpolation for algorithm
		up[k] = (1.0/(dx*dy))*	(u[ij]*((i+1)*dx - xp[k])*((j+1)*dy - yp[k]) + 
								u[ije]*(xp[k] - i*dx)*((j+1)*dy - yp[k]) + 
								u[ijn]*((i+1)*dx - xp[k])*(yp[k] - j*dy) +
								u[ijne]*(xp[k] - i*dx)*(yp[k] - j*dy));
		
		vp[k] = (1.0/(dx*dy))*	(v[ij]*((i+1)*dx - xp[k])*((j+1)*dy - yp[k]) + 
								v[ije]*(xp[k] - i*dx)*((j+1)*dy - yp[k]) + 
								v[ijn]*((i+1)*dx - xp[k])*(yp[k] - j*dy) +
								v[ijne]*(xp[k] - i*dx)*(yp[k] - j*dy));
		
	}
	
}