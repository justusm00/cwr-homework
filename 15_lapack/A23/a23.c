#include<tgmath.h>
#include<stdlib.h>
#include "my_numerics.h"
#include<stdio.h>


int main()
{
	FILE *matrix = fopen("23_matrix.dat","r"); //open matrix file
	FILE *vector = fopen("23_b.dat","r"); //open vector file
	FILE *mod1 = fopen("modulus1.dat","w");
	FILE *mod2 = fopen("modulus2.dat","w");
	double mod_jac = 0; //jacobi modulus
	double mod_gs = 0; //gauss seidel modulus
	
	int n = 10; //matrix dimension
	double A[n*n];
	double b[n];
	for(int i = 0; i < n*n; i++) fscanf(matrix, "%lf", &A[i]); //fill matrix
	for(int i = 0; i< n; i++) fscanf(vector, "%lf", &b[i]); //fill vector
	
	double x_jac[n]; //solution vector for jacobi method
	double x_gs[n]; //solution vector for gauss seidel method
	for(int i = 0; i < n; i++) 
	{
		x_jac[i] = 0; //initial guess, I guess
		x_gs[i] = 0;
	}
	int k_max = 40; //no of iterations
	
	//solution using jacobi
	lgs_jacobi_solve(n,A,b,x_jac,k_max);
	printf("Solution using Jacobi method: \n");
	for(int i = 0; i < n; i++) printf("%lf\n", x_jac[i]);
	
	

	
	//solution using gauss seidel
	lgs_gs_solve(n,A,b,x_gs,k_max);
	printf("Solution using gauss seidel method: \n");
	for(int i = 0; i < n; i++) printf("%lf\n", x_gs[i]);
	
	
	
	for(int i = 0; i < n; i++) 
	{
		x_jac[i] = 0; //initial guess, I guess
		x_gs[i] = 0;
	}
	
	
	//modulus and modulus of difference vectors
	
	
	//difference vectors xk - xk+1
	double diff_jac[n]; 
	double diff_gs[n];
	
	//kth iteration
	double xk_jac[n];
	double xk_gs[n];
	
	
	
	for(int i = 0; i < k_max; i++)
	{
	
		//modulus of xk
		for(int j = 0; j < n; j++)
		{
			mod_jac += x_jac[j] * x_jac[j]; //add to modulus
			mod_gs += x_gs[j] * x_gs[j];
			xk_jac[j] = x_jac[j]; //copy x to xk
			xk_gs[j] = x_gs[j];
			
		}
		mod_jac = sqrt(mod_jac);
		mod_gs = sqrt(mod_gs);
		fprintf(mod1,"%d\t%.10lf\t%.10lf\n",i,mod_jac,mod_gs); //print modulus
		lgs_jacobi_solve(n,A,b,x_jac,1); //do one iteration
		lgs_gs_solve(n,A,b,x_gs,1);
		mod_jac = 0;
		mod_gs = 0;
		
		//modulus of xk-xk+1
		//compute xk-xk+1
		for(int j = 0; j < n; j++)
		{
			diff_jac[j] = xk_jac[j] - x_jac[j];
			diff_gs[j] = xk_gs[j] - x_gs[j];
		
		}
		
		
		//compute modulus
		for(int j = 0; j < n; j++)
		{
			mod_jac += diff_jac[j] * diff_jac[j]; //add to modulus
			mod_gs += diff_gs[j] * diff_gs[j];
			
		}
		mod_jac = sqrt(mod_jac);
		mod_gs = sqrt(mod_gs);
		fprintf(mod2,"%d\t%.10lf\t%.10lf\n",i,mod_jac,mod_gs); //print modulus
		mod_jac = 0;
		mod_gs = 0;
	
	}
	for(int i = 0; i < n; i++) 
	{
		x_jac[i] = 0; //initial guess, I guess
		x_gs[i] = 0;
	}
	
	
	
	
	fclose(matrix);
	fclose(vector);
	fclose(mod1);
	fclose(mod2);
	return 0;



}
