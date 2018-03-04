#include "cp1headers.h"

void writeflowfield(int nx, int ny, float dx, float dy, float* u, float* v, float* p){
	
	float x,y;
	int i,j,ij;
	FILE *flowfield;
	
	flowfield = fopen("velocity.plt", "w");
	
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

void writeparticlepos(float* xp, float* yp, float* up, float* vp, int step, int np){
	
	int i;
	char filename[100];
	FILE *particlepos;
	
	sprintf(filename,"pos%d.plt",step);
	
	particlepos = fopen(filename,"w");
	
	fprintf(particlepos,"xp,yp,up,vp\n");
	
	for(i=0; i<np; i++){
		fprintf(particlepos,"%f,%f,%f,%f\n",xp[i],yp[i],up[i],vp[i]);
	}
	
	fclose(particlepos);
	
}