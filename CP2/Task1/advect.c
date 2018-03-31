#include "cp2task1headers.h"

void advectx(int nx, int ny, float dx, float dt, float* u, float* phi, float* phi_new, float epsilon){
	/*Advects in the x-direction using first order upwinding*/
	int i,j,ij,ije,ijw;
	
	for(i=1;i<nx+1;i++){
		for(j=1;j<ny+1;j++){
			ij = i +j*(nx+2);
			ije = (i+1)+ j*(nx+2);
			ijw = (i-1)+ j*(nx+2);
			
			
			//if(fabs(phi[ij])<=epsilon)	
				if(1)
			{
				if(u[ij]>=0.0){
					phi_new[ij] = phi[ij] - dt*u[ij]*(phi[ij]-phi[ijw])/dx;
				}
				else{
					phi_new[ij] = phi[ij] - dt*u[ij]*(phi[ije]-phi[ij])/dx;
				}
			}
			
		}
	}
}

void advecty(int nx, int ny, float dy, float dt, float* v, float* phi, float* phi_new, float epsilon){
	/*Advects in the y-direction using first order upwinding*/
	int i,j,ij,ijn,ijs;
	
	for(i=1;i<nx+1;i++){
		for(j=1;j<ny+1;j++){
			ij = i +j*(nx+2);
			ijn = i + (j+1)*(nx+2);
			ijs = i + (j-1)*(nx+2);
			
					
			//if(fabs(phi[ij])<=epsilon)	
				if(1)
			{
				if(v[ij]>=0.0){
					phi_new[ij] = phi_new[ij] - dt*v[ij]*(phi[ij]-phi[ijs])/dy;
				}
				else{
					phi_new[ij] = phi_new[ij] - dt*v[ij]*(phi[ijn]-phi[ij])/dy;
				}
			}
			
		}
	}
}