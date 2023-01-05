#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define PI 3.14159265358979323846

double integrand(double x);
double erf_midpoint(double x, double delta_x);
double erf_simpson(double x, double delta_x);


int main()
{
	//Berechnung mit der Simpson Regel im Interval [-2,2]
	FILE* myFile = fopen("data3.txt", "w");
	int N = 100;
	double x = -2.0;
	double i = 4.0 / N;
	double delta_x = 1e-4;
	for(int j = 0; j <= N; j++)
	{
			fprintf(myFile, "%.15lf \t %.15lf \n", x, erf_simpson(x, delta_x));
			x += i;
	}
	
	return 0;
}


double erf_midpoint(double x, double delta_x)
{
	int N = x / delta_x;
	double result = 0;
	double m = 0;
	for(int i = 1; i <= N; i++){
		m = (i * delta_x + (i + 1) * delta_x) / 2;
		result += integrand(m);
	}
	result *= delta_x;
	return result;
}


double erf_simpson(double x, double delta_x)
{
	int sign = 1;
	if(x < 0) sign = -1;
	int N = sign * x / delta_x;
	double result = 0;
	double m = 0;
	for(int i = 1; i <= N; i++){
		m = (i * delta_x + (i + 1) * delta_x) / 2;
		result += integrand(i * delta_x) + 4 * integrand(m) + integrand((i + 1) * delta_x);
	}
	result *= delta_x/6 * sign;
	return result;
}


double integrand(double x)
{
	double result = 2/sqrt(PI) * exp(-x * x);
	return result;
}
