#include "cp2task1headers.h"

void bcvelocity(int nx, int ny, float* u, float* v){
	
	/*Sets zero-gradient collocated BCs*/
	int i,j,ij,ijw,ije,ijn,ijs;
	
	for(i=0;i<nx+2;i++){
		//Bottom surface
		j=0;
		ij  = i +j*(nx+2); 
		ijn = i + (j+1)*(nx+2);
		
		u[ij] = u[ijn];
		v[ij] = v[ijn];
		
		//Top surface
		j = ny+1;
		ij  = i +j*(nx+2);
		ijs = i + (j-1)*(nx+2);
		
		u[ij] = u[ijs];
		v[ij] = v[ijs];
	}
	
	for(j=0;j<ny+2;j++){
		//Left surface
		i=0;
		ij  = i +j*(nx+2); 
		ije = (i+1)+ j*(nx+2);
		
		u[ij] = u[ije];
		v[ij] = v[ije];
		
		//Right surface
		i=nx+1;
		ij  = i +j*(nx+2); 
		ijw = (i-1)+ j*(nx+2);
		
		u[ij] = u[ijw];
		v[ij] = v[ijw];
		
	}
}

void bcphi(int nx, int ny, float* phi){
	/*Sets zero-gradient collocated BCs*/
	int i,j,ij,ijw,ije,ijn,ijs;
	
	for(i=0;i<nx+2;i++){
		//Bottom surface
		j=0;
		ij  = i +j*(nx+2); 
		ijn = i + (j+1)*(nx+2);
		
		phi[ij] = phi[ijn];
		
		//Top surface
		j = ny+1;
		ij  = i +j*(nx+2);
		ijs = i + (j-1)*(nx+2);
		
		phi[ij] = phi[ijs];
	}
	
	for(j=0;j<ny+2;j++){
		//Left surface
		i=0;
		ij  = i +j*(nx+2); 
		ije = (i+1)+ j*(nx+2);
		
		phi[ij] = phi[ije];
		
		//Right surface
		i=nx+1;
		ij  = i +j*(nx+2); 
		ijw = (i-1)+ j*(nx+2);
		
		phi[ij] = phi[ijw];
		
	}
	
}