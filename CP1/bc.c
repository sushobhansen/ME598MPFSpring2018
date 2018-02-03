#include "cp1headers.h"

void bc(int nx, int ny, float* u, float* v, float utop, bool intermediate){
	
	int i,j,ij,ijw,ije,ijn,ijs;
	
	for ( j = 1; j < ny + 1; j++){
      i = 0;
      ij  = i +j*(nx+2); ije = (i+1)+ j*(nx+2);
      u[ij] = -u[ije];
      v[ij] = -v[ije];

      i = nx+1;
      ij  = i +j*(nx+2); ijw = (i-1)+ j*(nx+2);
      u[ij] = -u[ijw];
      v[ij] = -v[ijw];
    }
	
    for ( i = 1; i < nx + 1; i++){
		j = 0;
		ij  = i +j*(nx+2); ijn = i + (j+1)*(nx+2);
		v[ij] = -v[ijn];
		u[ij] = -u[ijn];

		j = ny+1;
		ij  = i +j*(nx+2);ijs = i + (j-1)*(nx+2);

		if(intermediate){
		  v[ij] = -v[ijs];
		  u[ij] = -u[ijs];
		}
		else{
			u[ij] = 2*utop - u[ijs];
			v[ij] = -v[ijs];
		}
    }
}