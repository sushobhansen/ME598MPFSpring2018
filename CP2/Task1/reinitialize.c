#include "cp2task1headers.h"

void reinitialize(int nx, int ny, int niter, float dx, float dy, float dtau, float* phi, float* phi_new){
	int i,j,ij,ijn,ijs,ije,ijw,iter;
	float Dx,Dy;
	
	for(iter=0;iter<niter;iter++){
		printf("Reinitialization at iter = %d\n",iter+1);
		for(i=1;i<nx+1;i++){
			for(j=1;j<ny+1;j++){
				ij = i +j*(nx+2);
				ijn = i + (j+1)*(nx+2);
				ijs = i + (j-1)*(nx+2);
				ije = (i+1)+ j*(nx+2);
				ijw = (i-1)+ j*(nx+2);
				
				Dx = (phi[ije]-phi[ijw])/(2.0*dx);
				Dy = (phi[ijn]-phi[ijs])/(2.0*dy);
				
				phi_new[ij] = phi[ij] + dtau*sign(phi[ij])*(1.0-sqrt(pow(Dx,2.0)+pow(Dy,2.0)));
				swap(nx, ny, phi, phi_new); //Copy phi_new into phi, phi is lost
			}
		}
	}
}

float sign(float phi){
	return phi/fabs(phi);
}