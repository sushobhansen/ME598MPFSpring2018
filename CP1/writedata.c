#include "cp1headers.h"

void writedata(int nx, int ny, float dx, float dy, float* u, float* v, float* p){
	
	float x,y;
	int i,j,ij;
	FILE *flowfield;
	
	flowfield = fopen("driven_cavity_flow_c.plt", "w");
	
	fprintf(flowfield,"x,y,u,v,p\n");
	
	for ( j = 1; j < ny + 1; j++){
		for ( i = 1; i < nx + 1; i++){
			ij = i + j*(nx+2);
			x = (i - 1) * dx;
			y = (j - 1) * dy;
			fprintf(flowfield, "%f,%f,%f,%f,%f\n", x, y, u[ij], v[ij], p[ij]);
		}
	}
	
	fclose(flowfield);
	
}