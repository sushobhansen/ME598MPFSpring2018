#include "cp2task1headers.h"

void advectx(int nx, int ny, float dx, float dt, float* u, float* phi, float* phi_new, float epsilon){
	/*Advects in the x-direction using first order upwinding*/
	int i,j,ij,ije,ijw;
	
	for(i=1;i<nx+1;i++){
		for(j=1;j<ny+1;j++){
			ij = i +j*(nx+2);
			ije = (i+1)+ j*(nx+2);
			ijw = (i-1)+ j*(nx+2);
			
			if(fabs(phi[ij])<=epsilon)	
				//if(1)
			{
				//First order upwinding
				if(u[ij]>=0.0){
					phi_new[ij] = phi[ij] - dt*u[ij]*(phi[ij]-phi[ijw])/dx;
				}
				else{
					phi_new[ij] = phi[ij] - dt*u[ij]*(phi[ije]-phi[ij])/dx;
				}
				
				/* //QUICK
				if(u[ij]>=0.0){
					phi_new[ij] = phi[ij] - dt*u[ij]*((3.0*phi[ije]+6.0*phi[ij]-phi[ijw])-(3.0*phi[ij]+6.0*phi[ijw]-phi[ijww]))/(8.0*dx);
				}
				else{
					phi_new[ij] = phi[ij] - dt*u[ij]*((3.0*phi[ij]+6.0*phi[ije]-phi[ijee])-(3.0*phi[ijw]+6.0*phi[ij]-phi[ije]))/(8.0*dx);
				}  */
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
			
					
			if(fabs(phi[ij])<=epsilon && j!=1 && j!=ny)	
				//if(1)
			{
				//First order upwinding
				if(v[ij]>=0.0){
					phi_new[ij] = phi_new[ij] - dt*v[ij]*(phi[ij]-phi[ijs])/dy;
				}
				else{
					phi_new[ij] = phi_new[ij] - dt*v[ij]*(phi[ijn]-phi[ij])/dy;
				}
				
				//QUICK				
				/* if(v[ij]>=0.0){
					phi_new[ij] = phi_new[ij] - dt*v[ij]*((3.0*phi[ijn]+6.0*phi[ij]-phi[ijs])-(3.0*phi[ij]+6.0*phi[ijs]-phi[ijss]))/(8.0*dy);
				}
				else{
					phi_new[ij] = phi_new[ij] - dt*v[ij]*((3.0*phi[ij]+6.0*phi[ijn]-phi[ijnn])-(3.0*phi[ijs]+6.0*phi[ij]-phi[ijn]))/(8.0*dy);
				} */
			}
			
		}
	}
}