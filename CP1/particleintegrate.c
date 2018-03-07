#include "cp1headers.h"

/*
Integrates velocity using 4th Order Runge-Kutta scheme
*/
void integratevelocity(float* uf, float* up_old, float* up_new, float dtp, float U, float St, float L, int np){
	
	int k; 
	float K1, K2, K3, K4;
	
	for (k=0; k<np; k++){
		
		K1 = drag(uf[k],up_old[k], U, St, L);
		K2 = drag(uf[k],up_old[k] + dtp*K1*0.5, U, St, L);
		K4 = drag(uf[k],up_old[k] + dtp*K3, U, St, L);
		K3 = drag(uf[k],up_old[k] + dtp*K2*0.5, U, St, L);
		
		up_new[k] = up_old[k] + dtp*(K1+2.0*K2+2.0*K3+K4)/6.0;
		
	}
	
}

/*
Integrates xp using the trapezoidal rule
*/
void integrateposition(float* xp, float* up_old, float* up_new, float dtp, int np){
	int k;
	
	for (k=0; k<np; k++){
		xp[k] = xp[k] + 0.5*dtp*(up_old[k]+up_new[k]);
		if(xp[k]>1.0){xp[k]=1.0;} //If particle leaves domain, put it back in the domain
	}
	
}

/*
Evaluates drag force per unit mass
*/
float drag(float uf, float up, float U, float St, float L){
	//Drag force per unit mass for Stokes flow, dup/dt = (U/(St*L))*(uf-up)
	return (U/(St*L))*(uf-up);
}

/*
Puts the values of up_new into up_old. Values of up_old are lost.
*/
void swap(float* up_old, float* up_new, int np){
	int k;
	
	for (k=0; k<np; k++){
		up_old[k] = up_new[k];		
	}
	
}