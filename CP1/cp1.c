#include "cp1headers.h"

#define nx 512
#define ny 512

int main (int argc, char *argv[]){
	
	float *u, *uh, *v, *vh, *p, *an, *as, *ae, *aw, *ap, *s, *ttemp, *ph;
	float x1, y1, dx, dy, dt, utop, Re, amu, omega;
	int nt, niter, it, i, j, size, ij, ije, ijw, ijs, ijn;
	clock_t start, end; 
	
	size = (nx+2)*(ny+2)*sizeof(float);
	
	//Allocate memory to each array
	u  = (float*) malloc (size);
	uh = (float*) malloc (size);
	v  = (float*) malloc (size);
	vh = (float*) malloc (size);
	p  = (float*) malloc (size);
	ph = (float*) malloc (size);
	an = (float*) malloc (size);
	as = (float*) malloc (size);
	ae = (float*) malloc (size);
	aw = (float*) malloc (size);
	ap = (float*) malloc (size);
	s  = (float*) malloc (size);
	
	//Define dimensions of the geometry
	x1 = 1.0;                 // Length of domain along x - axis
	y1 = 1.0;                 // Length of domain along y - axis
	
	//Flow variables
	dx = x1/nx;
	dy = y1/ny;
	dt = 0.000075;                // Time step
	utop = 1.0;               	// Velocity of the Lid
	Re = 100.0;               // Reynolds number of the flow
	amu = x1 * utop / Re;     // We are calculating the viscosity from the Reynolds number of the flow
	nt = 100;                // Maximum number of iterations
	niter = 200;               // Number of iteration for the Pressure Poisson solver
	omega = 1.0;              // Optimum value of the omega for SOR method
	
	//Initialize all the variables to zero everywhere
	initialize(nx, ny, u, v, p);
	
	//Apply the boundary condition on u on the top wall
	bc(nx, ny, u, v, utop);
	
	//Evaluate the coefficients of the Pressure Poisson Equation (PPE)
	ppecoeffs(nx, ny, aw, ae, an, as, ap, dx, dy); 
	
	//Start time march
	start = clock();
	
	for ( it = 0; it < nt; it ++){
		printf("time step = %d\n", it+1);
		
		//Solve for uhat and vhat, apply BCs and evaluate source term for PPE
		momentum(nx, ny, u, v, uh, vh, s, amu, dx, dy, dt);
		bc(nx, ny, uh, vh, utop);
		ppesource(nx, ny, uh, vh, s, dx, dy, dt);
		
		//Run SOR iterations to solve for p
		sor(nx, ny, niter, aw, ae, an, as, ap, p, s, ph, ttemp, omega);
		
		//Correct momentum for p and apply BCs
		momentumcorr(nx, ny, u, v, p, uh, vh, dx, dy, dt);
		bc(nx, ny, u, v, utop);
	}

	end = clock();
	printf("Time taken by CPU = %10.8f sec\n", ((float) (end - start)) / CLOCKS_PER_SEC);
	
	//Write data file for plotting
	writedata(nx, ny, dx, dy, u, v, p);
	
	return 0;
}
