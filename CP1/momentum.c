#include "cp1headers.h"

/*=================momentum.c=================
 *Marches u and v for diffusion and convection (not PG), applies BCs, and evaluates source terms for PPE
*/

void momentum(int nx, int ny, float* u, float* v, float* uh, float* vh, float* s, float amu, float dx, float dy, float dt){
	
	float conv_x, conv_y, diff_x, diff_y;
	int i,j,ij,ijw,ije,ijn,ijs;
	
	for ( i = 1; i < nx + 1; i++){
		for (j = 1; j < ny + 1; j++){
				
			ij  = i +j*(nx+2);
			ije = (i+1)+ j*(nx+2);
			ijw = (i-1)+ j*(nx+2);
			ijs = i + (j-1)*(nx+2);
			ijn = i + (j+1)*(nx+2);
			
			//Solve for x-momentum
			conv_x = u[ij] * (u[ije] - u[ijw])/ (2.0 * dx);
			conv_y = v[ij] * (u[ijn] - u[ijs])/ (2.0 * dy);
			diff_x = amu * (u[ije] - 2.0*u[ij] + u[ijw])/(dx * dx);
			diff_y = amu * (u[ijn] - 2.0*u[ij] + u[ijs])/(dy * dy);
			uh[ij] = u[ij] + dt * (-conv_x -conv_y + diff_x + diff_y);
			
			//Solve for y-momentum
			conv_x = u[ij] * (v[ije] - v[ijw])/ (2.0 * dx);
			conv_y = v[ij] * (v[ijn] - v[ijs])/ (2.0 * dy);
			diff_x = amu * (v[ije] - 2.0*v[ij] + v[ijw])/(dx * dx);
			diff_y = amu * (v[ijn] - 2.0*v[ij] + v[ijs])/(dy * dy);
			vh[ij] = v[ij] + dt * (-conv_x -conv_y + diff_x + diff_y);
		}
	}	
	
}