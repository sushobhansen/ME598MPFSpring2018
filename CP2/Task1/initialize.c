#include "cp2task1headers.h"

void initialize(int nx, int ny, float dx, float dy, float xcenter0, float ycenter0, float radius, float velocity, float* u, float* v, float* phi){
	
	int i,j,ij;
	float x,y;
	
	//Set values inside the domain
	for ( i = 1; i < nx + 1; i++){
		for ( j = 1; j < ny + 1; j++){
		  ij = i +j*(nx+2);
		  
		  //Set speed
		  u[ij] = velocity;
		  v[ij] = velocity;
		  
		  //Set signed distance - positive inside, negative outside
		  x = (float)(i-1)*dx;
		  y = (float)(j-1)*dy;
		  phi[ij] = radius-sqrt(pow(xcenter0-x,2)+pow(ycenter0-y,2));
		}
	}
	
}
