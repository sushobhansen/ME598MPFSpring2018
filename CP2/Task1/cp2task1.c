#include "cp2task1headers.h"

#define nx 321
#define ny 321

int main (int argc, char *argv[]){
	
	float *u, *v, *phi, *phi_new, *phi0;
	float xl, yl, dx, dy, dt, dtau, velocity, epsilon;
	float xcenter0, ycenter0, radius, xcenterf, ycenterf;
	int nt, niter, it, i, j, size, ij, ije, ijw, ijs, ijn;
	clock_t start, end; 
	
	size = (nx+2)*(ny+2)*sizeof(float);
	
	//Allocate memory to each array
	u  = (float*) malloc (size);
	v  = (float*) malloc (size);
	phi = (float*) malloc (size);
	phi_new = (float*) malloc (size);
	phi0 = (float*) malloc (size);
	
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
	dt = 0.5/(velocity/dx + velocity/dy); //Courant number of 0.5 for main time step
	dtau = 0.1/(velocity/dx + velocity/dy); //Courant number of 0.1 for reinitialization
	nt = (int)ceil((fabs(xcenterf-xcenter0)/velocity)/dt);                // Number of iterations to bring circle to final position
	niter = 10;               // Number of iteration for reinitialization
	epsilon = 30.0*dx; //bandwidth
	
	//Initialize and set BCs (zero-gradient) for velocity and phi
	initialize(nx, ny, dx, dy, xcenter0, ycenter0, radius, velocity, u, v, phi);
	bcvelocity(nx, ny, u, v);
	bcphi(nx,ny,phi);
	swap(nx, ny, phi_new, phi); //Copy phi into phi_new
	
	//Write initial data
	writeflowfield(nx, ny, dx, dy, u, v);
	writephi(nx, ny, dx, dy, phi, 0);
	
	//Begin iterating
	printf("Beginning advection, dt = %f, nt = %d, dx = %f, dy = %f\n",dt,nt,dx,dy);
	start = clock();
	
	for (it = 0; it < nt; it++){
		printf("Advection time step = %d\n",it+1);
		advectx(nx, ny, dx, dt, u, phi, phi_new, epsilon); //advects in x, phi* stored in phi_new
		bcphi(nx, ny, phi_new); //set bc for phi_new
		advecty(nx, ny, dy, dt, v, phi, phi_new, epsilon); //advects in y, phi_new stores final value
		bcphi(nx, ny, phi_new); //set bc for phi_new
		swap(nx, ny, phi, phi_new); //Copy phi_new into phi, phi is advected field
				
		printf("Reinitializing...\n");
		reinitialize(nx, ny, niter, dx, dy, dtau, phi, phi_new, phi0, epsilon); //reinitialize, phi_new is reinitialized field
		swap(nx, ny, phi, phi_new); //phi is reinitialized field
		bcphi(nx,ny,phi);
		
		if((it+1)%80==0){
		writephi(nx, ny, dx, dy, phi, it+1); 
		printf("Writing to file...\n");
		}
	}
	
	end = clock();
	
	printf("Time taken by CPU for solution = %10.8f sec\n", ((float) (end - start)) / CLOCKS_PER_SEC);
	printf("Simulation time = %f, time step = %f\n",dt*(float)nt,dt);
	
	
	return 0;
}
