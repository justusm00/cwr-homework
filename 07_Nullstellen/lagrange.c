#include<stdio.h>
#include<stdlib.h>
#include<tgmath.h>
#include "my_numerics.h"


const double mu_earth = 3.986e14;
const double mu_sun = 1.327e20;
const double omega = 1.991e-7;
const double R = 1.496e11;

double acceleration(double r);



int main(){
	double r = find_root(acceleration, R + 100, 1000, 10, 1000);
	printf("Der Lagrangepunkt L2 ist ca %g m weit von der Erde entfernt\n", r - R);
}


double acceleration(const double r)
{
	double a_earth = -mu_earth / ((r-R) * (r - R));
	double a_sun = -mu_sun / (r * r);
	double a_centrifugal = omega * omega * r;
	return a_earth + a_sun + a_centrifugal;
}



