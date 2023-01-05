#include<stdio.h>
#include<stdlib.h>
#include<tgmath.h>
#include "my_numerics.h"


double f(double x);

int main()
{
	int N = 100;
	double delta = 1.0;
	double x = 1.0;
	double res = 1.0; //analytische Lösung für x = 1
	FILE* myFile = fopen("derive.txt", "w");
	for(int i = 0; i < N; i++)
	{
		fprintf(myFile, "%.16lf \t %.16lf \n", delta , fabsf(diff(x, delta, f) - res));
		delta *= 0.7;
	}
	return 0;


}



double f(double x){
	return x * x * (x - 1);
}
