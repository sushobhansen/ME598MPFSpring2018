#include "cp2task1headers.h"

void advect(int nx, int ny, float dx, float dy, float dt, float* u, float* v, float* phi, float* phi_new){
	
	int i,j,ij,ije,ijn,ijw,ijs;
	float k, F;
	float deltap, deltam, Dpx, Dmx, Dpy, Dmy;
	
	k = 3.0*dx; //bandwidth
	
	for(i=1;i<nx+1;i++){
		for(j=1;j<ny+1;j++){
			ij = i +j*(nx+2);
			ijn = i + (j+1)*(nx+2);
			ijs = i + (j-1)*(nx+2);
			ije = (i+1)+ j*(nx+2);
			ijw = (i-1)+ j*(nx+2);
			
			if(fabs(phi[ij])<=k){
				F = sqrt(pow(u[ij],2.0)+pow(v[ij],2.0));
				
				//Calculate delta's
				Dpx = (phi[ije]-phi[ij])/dx;
				Dmx = (phi[ij]-phi[ijw])/dx;
				Dpy = (phi[ijn]-phi[ij])/dy;
				Dmy = (phi[ij]-phi[ijs])/dy;
				
				deltap = sqrt(pow(fmax(Dmx,0.0),2.0)+pow(fmin(Dpx,0.0),2.0)+pow(fmax(Dmy,0.0),2.0)+pow(fmin(Dpy,0.0),2.0));
				deltam = sqrt(pow(fmax(Dpx,0.0),2.0)+pow(fmin(Dmx,0.0),2.0)+pow(fmax(Dpy,0.0),2.0)+pow(fmin(Dmy,0.0),2.0));
				
				//Integrate
				phi_new[ij] = phi[ij] - dt*(fmax(F,0.0)*deltap+fmin(F,0.0)*deltam);
			}
		}
	}
	
}