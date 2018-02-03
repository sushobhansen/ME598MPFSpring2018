#include "cp1headers.h"

void ppecoeffs(int nx, int ny, float* aw, float* ae, float* an, float* as, float* ap, float dx, float dy){
	int i,j,ij,ijw,ije,ijn,ijs;
	
	for ( i = 1; i < nx+1; i++){
		for ( j = 1; j < ny + 1; j++){
		  ij  = i +j*(nx+2);
		  ije = (i+1)+ j*(nx+2);
		  ijw = (i-1)+ j*(nx+2);
		  ijs = i + (j-1)*(nx+2);
		  ijn = i + (j+1)*(nx+2);
		  
		  
		  aw[ij] = 1.0/(dx * dx);
		  ae[ij] = 1.0/(dx * dx);
		  as[ij] = 1.0/(dy * dy);
		  an[ij] = 1.0/(dy * dy);
		  
		  //Remove links to ghost points to enforce dp/dn = 0
		  if ( i == 1 ) aw[ij] = 0.0;
		  if ( i == nx )ae[ij] = 0.0;
		  if ( j == 1 ) as[ij] = 0.0;
		  if ( j == ny )an[ij] = 0.0;
		  
		  ap[ij] = as[ij] + an[ij] + aw[ij] + ae[ij];
		}
	} 
	
}