#include "cp2task1headers.h"

#define nx 41
#define ny 41

int main (int argc, char *argv[]){
	
	float *u, *v, *phi, *phi_new;
	float xl, yl, dx, dy, dt, velocity;
	float xcenter0, ycenter0, radius, xcenterf, ycenterf;
	int nt, niter, it, i, j, size, ij, ije, ijw, ijs, ijn;
	clock_t start, end; 
	
	size = (nx+2)*(ny+2)*sizeof(float);
	
	//Allocate memory to each array
	u  = (float*) malloc (size);
	v  = (float*) malloc (size);
	phi = (float*) malloc (size);
	phi_new = (float*) malloc (size);
	
	//Define dimensions of the geometry
	xl = 8.0;                 // Length of domain along x - axis
	yl = 8.0;                 // Length of domain along y - axis
	xcenter0 = 2.0;			//Initial x-center of the circle
	ycenter0 = 2.0;			//Initial y-center of the circle
	radius = 0.5;			//Radius of the circle
	xcenterf = 6.0;			//Final x-center of the circle
	ycenterf = 6.0;			//Final y-center of the circle
	
	//Flow variables
	dx = xl/(float)(nx-1);
	dy = yl/(float)(ny-1);
	velocity = 1.0;               	// Fixed velocity components
	dt = 0.5/(velocity/dx + velocity/dy); //Courant number of 0.5
	nt = (int)ceil((fabs(xcenterf-xcenter0)/velocity)/dt);                // Number of iterations to bring circle to final position
	niter = 5;               // Number of iteration for reinitialization
	
	//Initialize and set BCs (zero-gradient) for velocity and phi
	initialize(nx, ny, dx, dy, xcenter0, ycenter0, radius, velocity, u, v, phi);
	bc(nx, ny, u, v, phi);
	swap(nx, ny, phi_new, phi); //Copy phi into phi_new
	
	//Write initial data
	writeflowfield(nx, ny, dx, dy, u, v);
	writephi(nx, ny, dx, dy, phi, 0);
	
	//Begin iterating
	printf("Beginning advection, dt = %f, nt = %d, dx = %f, dy = %f\n",dt,nt,dx,dy);
	start = clock();
	
	for (it = 0; it < nt; it++){
		printf("Advection time step = %d\n",it+1);
		advect(nx, ny, dx, dy, dt, u, v, phi, phi_new); //phi_new is the advected field
		swap(nx, ny, phi, phi_new); //Copy phi_new into phi, phi is lost
				
		printf("Reinitializing...\n");
		reinitialize(nx, ny, niter, dx, dy, dt*2.0, phi, phi_new);
		swap(nx, ny, phi, phi_new);
		bc(nx, ny, u, v, phi);
		writephi(nx, ny, dx, dy, phi, it+1); printf("Writing to file...\n");
	}
	
	end = clock();
	
	printf("Time taken by CPU for solution = %10.8f sec\n", ((float) (end - start)) / CLOCKS_PER_SEC);
	printf("Simulation time = %f, time step = %f\n",dt*(float)nt,dt);
	
	
	return 0;
}
