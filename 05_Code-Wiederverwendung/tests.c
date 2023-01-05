#include<stdio.h>
#include<stdlib.h>
#include "my_numerics.h"

int main(){

	//Test für solve_quadratic

	tuple x; 
	tuple x_ana; 
	tuple error; 
	x_ana.x1 = 1e-8;
	x_ana.x2 = 1e8;	
			
	double a = 1; 
	double b = (-1)*(1e16 + 1) / (1e8);
	double c = 1;
	x = solve_quadratic(a,b,c);
	error.x1 = (x_ana.x1 - x.x1) / x_ana.x1; 
	error.x2 = (x_ana.x2 - x.x2) / x_ana.x2;
	
	printf("Die Lösungen sind %.15lf und %.15lf \n", x.x1, x.x2);	
	printf("Die relativen Fehler sind %.15lf und %.15lf \n", error.x1, error.x2);
	
	
	//Test für mn_erf_simpson und mn_erf_mitpoint
	double delta_x = 1e-3;
	printf("erf(1.5) ist mit der Simpson-Regel ungefähr %.15lf \n", mn_erf_simpson(1.5,delta_x));
	printf("erf(1.5) ist mit der Midpoint-Regel ungefähr %.15lf \n", mn_erf_midpoint(1.5,delta_x));
	return 0;


}
