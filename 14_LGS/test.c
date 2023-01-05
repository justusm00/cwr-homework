#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include"my_numerics.h"
#include"20_matrix.h"

int main()
{	
	int n = 10;
	double B[n*n];
	double B_inverse[n*n];
	double x[n];
	double b[] = {1,0,0,0,0,0,0,0,0,0};
	for(int i = 0; i < n*n; i++)
	{
		B[i] = A[i];
	}
	
	
	
	mat_inverse(n, B, B_inverse);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			printf("%g\t", B_inverse[i * n + j]);
		}
		printf("\n");
	}
	
	
	return 0;
	
	
	
}
