#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<my_numerics.h>
#define PI 3.14159265358979323846

double lam_anal(int l, double w, double t, int n)
{
	return w + 2 * t * cos(l * PI / (n + 1));
}

double p(int n, double a[n], double b[n-1], double lambda, int grade)
{
	double p;
	double w = a[0];
	double t = -b[0];
	double p0 = 1; //starting values for recursive formula
	double p1 = w - lambda;
	for(int i = 0; i < n - 1; i++)
	{
		p = (w - lambda) * p1 - t * t * p0; //do recursion
		p0 = p1;
		p1 = p;
	}
	return p;

}


//need to make a function that depends only on lambda so we can use find_root -.-
double p_ev(double lambda)
{
	int n = 10;
	double t = 0.5; 
	double w = 1;
	double a[n]; //diagonal term
	double b[n-1]; //off-diagonal term
	double res;
	for(int i = 0; i < n; i++)
	{
		a[i] = w;	//fill diagonal
	}
	for(int i = 0; i < n-1; i++)
	{
		b[i] = -t;	//fill off diagonal
	}
	res = p(n,a,b,lambda,n);
	return res;
}




int main()
{
	FILE* posfile = fopen("pos.dat","w"); //open output file
	FILE* eigenvals = fopen("eigenvals.dat","w"); //open file to print eigenvalues to
	int n = 10;
	double t = 0.5; 
	double w = 1;
	
	double lambda[n]; //store results here
	
	//set boundaries according to analytical solution
	double min = w - 2* t;
	double max = w + 2* t;
	double step = 1e-2;
	double x = min;
	double res;
	while(x <= max)
	{
		res = p_ev(x);
		fprintf(posfile,"%g \t %g\n", x,res);
		x += step;
	}
	
	//compute eigenvalues and save them to lambda
	double x0[] = {0.01, 0.1, 0.3, 0.55, 0.85, 1.1, 1.4, 1.6, 1.8, 1.92}; //initial "guesses"
	for(int i = 0; i < n; i++)
	{
		lambda[i] = find_root(p_ev,x0[i] , 1e-4, 1e-4,1000);
		fprintf(eigenvals,"%d\t%.10g\t %.10g\n",n-i,lambda[i],lam_anal(n-i,w,t,n));
		printf("%g \n", lambda[i]);
	}
	
	
	
	
	return 0;


}
