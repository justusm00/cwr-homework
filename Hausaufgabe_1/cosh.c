#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define PI 3.14159265358979323846

double integrand(double x);
double cosh_midpoint(double x, double delta_x);
double cosh_simpson(double x, double delta_x);


int main()
{

	FILE* myFile = fopen("cosh.txt", "w");
	double value = sinh(1);
	int N = 21;
	double x = 1.0;
	double delta_x = 0.1;
	for(int j = 0; j < N; j++)
	{
			fprintf(myFile, "%.10lf \t %.20lf \t %.20lf \n", delta_x, fabsf(value - cosh_simpson(x, delta_x)), fabsf(value - cosh_midpoint(x, delta_x)));
			delta_x /= 2;
	}
	fclose(myFile);
	
}


double cosh_midpoint(double x, double delta_x)
{
	int N = x / delta_x;
	double result = 0;
	double m = 0;
	for(int i = 0; i < N; i++){
		m = (i * delta_x + (i + 1) * delta_x) / 2;
		result += integrand(m);
	}
	result *= delta_x;
	return result;
}


double cosh_simpson(double x, double delta_x)
{
	int sign = 1;
	if(x < 0) sign = -1;
	int N = sign * x / delta_x;
	double result = 0;
	double m = 0;
	for(int i = 0; i < N; i++){
		m = (i * delta_x + (i + 1) * delta_x) / 2;
		result += integrand(i * delta_x) + 4 * integrand(m) + integrand((i + 1) * delta_x);
	}
	result *= delta_x/6 * sign;
	return result;
}


double integrand(double x)
{
	double result = cosh(x);
	return result;
}
