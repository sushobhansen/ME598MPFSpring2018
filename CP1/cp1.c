#include "cp1headers.h"

#define nx 41
#define ny 41

int main (int argc, char *argv[]){
	
	float *u, *uh, *v, *vh, *p, *an, *as, *ae, *aw, *ap, *s, *ttemp, *ph;
	float xl, yl, dx, dy, dt, utop, Re, amu, omega;
	int nt, niter, it, i, j, size, ij, ije, ijw, ijs, ijn;
	clock_t start, end; 
	
	//Variables for particle tracking
	float *up_old, *vp_old, *up_new, *vp_new, *uf, *vf, *xp, *yp;
	float St, dtp;
	int np, ntp, itp, nplot;
	
	size = (nx+2)*(ny+2)*sizeof(float);
	np = 50; //No of particles
	
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
	
	up_old = (float*) malloc (np*sizeof(float));
	vp_old = (float*) malloc (np*sizeof(float));
	up_new = (float*) malloc (np*sizeof(float));
	vp_new = (float*) malloc (np*sizeof(float));
	uf = (float*) malloc (np*sizeof(float));
	vf = (float*) malloc (np*sizeof(float));
	xp = (float*) malloc (np*sizeof(float));
	yp = (float*) malloc (np*sizeof(float));
	
	//Define dimensions of the geometry
	xl = 1.0;                 // Length of domain along x - axis
	yl = 1.0;                 // Length of domain along y - axis
	
	//Flow variables
	dx = xl/(float)(nx-1);
	dy = yl/(float)(ny-1);
	utop = 1.0;               	// Velocity of the Lid
	dt = 0.5*dx/utop;
	Re = 100.0;               // Reynolds number of the flow
	amu = xl * utop / Re;     // We are calculating the viscosity from the Reynolds number of the flow
	nt = 500;                // Maximum number of iterations
	niter = 200;               // Number of iteration for the Pressure Poisson solver
	omega = 0.8;              // Optimum value of the omega for SOR method
	
	
	
	//Initialize all the variables to zero everywhere
	initialize(nx, ny, np, u, v, p, up_old, vp_old);
	
	printf("****Solving for Velocity and Pressure Fields****\n");
	
	//Apply the boundary condition on u on the top wall
	bc(nx, ny, u, v, utop);
	
	//Evaluate the coefficients of the Pressure Poisson Equation (PPE)
	ppecoeffs(nx, ny, aw, ae, an, as, ap, dx, dy); 
	
	//Start time march
	start = clock();
	
	for ( it = 0; it < nt; it ++){
		printf("fluid time step = %d\n", it+1);
		
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
	printf("Time taken by CPU for fluid solution = %10.8f sec\n", ((float) (end - start)) / CLOCKS_PER_SEC);
	printf("Simulation time = %f, time step = %f\n",dt*(float)nt,dt);
	
	//Write velocity field file for plotting
	writeflowfield(nx, ny, dx, dy, u, v, p);
	
	printf("****Solving for Particle Positions****\n");
	//Parameters for particle tracking
	St = 1.0; //Stokes numbers
	dtp = dt; //Set time step for particles equal to time step for fluid (guess)
	ntp = 2000; //Set simulation time for particles equal to that for fluid (guess)
	nplot = 400; //Number of time steps to plot at
	
	//Randomly initialize np particles in [0.25,0.75]^2
	for (i=0; i<np; i++){
		xp[i] = 0.25 + (0.75-0.25)*(float)rand()/(float)RAND_MAX;
		yp[i] = 0.25 + (0.75-0.25)*(float)rand()/(float)RAND_MAX;		
	}
	
	//Interpolate fluid velocity to initial particle positions
	bilinterp(u, v, xp, yp, uf, vf, np, nx, ny, dx, dy);
	
	//Write initial particle position for plotting
	writeparticlepos(xp, yp, up_old, vp_old, 0, np);
	
	//Integrate from 0 to ntp
	itp = 0;
	start = clock();
	
	while(itp<ntp){
		printf("particle time step = %d\n",itp+1);
		
		//Integrate u and x
		integratevelocity(uf, up_old, up_new, dtp, utop, St, xl, np);
		integrateposition(xp, up_old, up_new, dtp, np);
		
		//Integrate v and y
		integratevelocity(vf, vp_old, vp_new, dtp, utop, St, xl, np);
		integrateposition(yp, vp_old, vp_new, dtp, np);
		
		//Write velocity at certain times
		if((itp+1)%nplot == 0){
			writeparticlepos(xp, yp, up_new, vp_new, itp+1, np);			
		}
		
		//Put new values into old values for next time step
		swap(up_old, up_new, np);
		swap(vp_old, vp_new, np);
		
		//Interpolate fluid velocity to new particle positions
		bilinterp(u, v, xp, yp, uf, vf, np, nx, ny, dx, dy);
		
		//Next time step
		itp++;
	}
	
	end = clock();
	
	printf("Time taken by CPU for particle tracking = %10.8f sec\n", ((float) (end - start)) / CLOCKS_PER_SEC);
	printf("Simulation time = %f, time step = %f\n",dtp*(float)ntp,dtp);
	
	
	
	return 0;
}
