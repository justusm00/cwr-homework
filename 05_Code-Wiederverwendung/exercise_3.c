#include<stdio.h>
#include<stdlib.h>
#include "my_numerics.h"
#define PI 3.14159265358979323846


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
			fprintf(myFile, "%.15lf \t %.15lf \n", x, mn_erf_simpson(x, delta_x));
			x += i;
	}
	
	return 0;
}



