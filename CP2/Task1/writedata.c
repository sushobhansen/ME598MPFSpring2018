#include "cp2task1headers.h"

void writeflowfield(int nx, int ny, float dx, float dy, float* u, float* v){
	
	float x,y;
	int i,j,ij;
	FILE *flowfield;
	
	flowfield = fopen("velocity.plt", "w");
	
	fprintf(flowfield,"x,y,u,v\n");
	
	for ( j = 1; j < ny + 1; j++){
		for ( i = 1; i < nx + 1; i++){
			ij = i + j*(nx+2);
			x = (i - 1) * dx;
			y = (j - 1) * dy;
			fprintf(flowfield, "%f,%f,%f,%f\n", x, y, u[ij], v[ij]);
		}
	}
	
	fclose(flowfield);
	
}

void writephi(int nx, int ny, float dx, float dy, float* phi, int step){
	
	int i,j,ij;
	float x,y;
	char filename[100];
	FILE *phifield;
	
	sprintf(filename,"phi%d.plt",step);
	
	phifield = fopen(filename,"w");
	
	fprintf(phifield,"x,y,phi\n");
	
	for ( j = 1; j < ny + 1; j++){
		for ( i = 1; i < nx + 1; i++){
			ij = i + j*(nx+2);
			x = (i - 1) * dx;
			y = (j - 1) * dy;
			fprintf(phifield, "%f,%f,%f\n", x, y, phi[ij]);
		}
	}
	
	fclose(phifield);
	
}

void swap(int nx, int ny, float* phi, float* phi_new){
	/*Swaps phi_new into phi and phi is lost, equivalent of phi = phi_new*/
	
	int i,j,ij;
	
	for(i = 0;i<nx+2;i++){
		for(j=0;j<ny+2;j++){
			ij = i +j*(nx+2);
			phi[ij] = phi_new[ij];			
		}
	}
	
}