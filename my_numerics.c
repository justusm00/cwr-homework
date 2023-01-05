#include<tgmath.h>
#include<stdlib.h>
#include "my_numerics.h"
#include<stdio.h>
#define PI 3.14159265358979323846

static double gaussian(double y)
{
return exp(-y*y);
}




double mn_integrate(double left, double right, int N, double integrand(double)){
	double width = right - left;
	double dx = width / N;
	double result = 0.0;
	for(int i = 0; i < N; i++){
		result += integrand(i * dx) * dx;
	}
	return result;
	
}





double mn_erf_simpson(double x, double delta_x)
{
	int sign = 1;
	if(x < 0) sign = -1;
	int N = sign * x / delta_x;
	double result = 0;
	double m = 0;
	for(int i = 1; i <= N; i++){
		m = (i * delta_x + (i + 1) * delta_x) / 2;
		result += gaussian(i * delta_x) + 4 * gaussian(m) + gaussian((i + 1) * delta_x);
	}
	result *= delta_x/6 * sign;
	result *= 2/sqrt(PI);
	return result;
}






double mn_erf_midpoint(double x, double delta_x)
{
	int N = x / delta_x;
	double result = 0;
	double m = 0;
	for(int i = 1; i <= N; i++){
		m = (i * delta_x + (i + 1) * delta_x) / 2;
		result += gaussian(m);
	}
	result *= delta_x;
	result *= 2/sqrt(PI);
	return result;
}






tuple solve_quadratic(double a, double b , double c)
{	
	tuple x;
	double d = sqrt(b * b - 4 * a * c);	
	double test1;
	double test2;
	
	//Vermeiden, dass Zähler bzw Nenner sehr klein werden
	if(fabsf(-b + d) > 1e-7)
	{
		x.x1 = (-b + d) / (2 * a);
		x.x2 = (2*c) / (-b + d);
	}
	else
	{
		x.x1 = (2*c) / (-b - d);
		x.x2 = (-b - d) / (2 * a);
	}
	test1 = x.x1;
	test2 = x.x2;
	
	//Ausgabe in aufsteigender Reihenfolge
	if(test1 > test2) 
	{
		x.x1 = test2; 
		x.x2 = test1;
	}
	return x;

}







double diff(double x, double delta, double func(double))
{
	return (func(x + delta) - func(x - delta)) / (2 * delta);
}



double find_root(double func(double), double x0, double delta, double rel_tol, int max_iter)
{	
	int i = 0;
	double x = x0;
	double newx = func(x0) / diff(x0,delta,func);
	while((i < max_iter) && (fabsf(newx - x) > rel_tol)){
		newx = x;
		x = x - func(x) / diff(x,delta,func);
	}
	return x;
		
}





void euler_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params)
{
	double *f = malloc((int)(sizeof(double) * dimension));
	func(t,y,f,params);
	for(int i = 0; i < dimension; i++)
	{
		y[i] += delta_t * (*(f + i));
	}
	free(f);
	


}



void rk2_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params)
{
	
	double *support = malloc((int)(sizeof(double) * dimension));
	double *k1 = malloc((int)(sizeof(double) * dimension));
	double *k2 = malloc((int)(sizeof(double) * dimension));
	func(t,y,k1,params);   //calculate k1
	for(int i = 0; i < dimension; i++)
	{
		k1[i] *= delta_t; //multiply each element by delta_t
	}
	for(int i = 0; i < dimension; i++)
	{
		support[i] = y[i] + 0.5 * k1[i]; //fill support array
	}
	func(t + 0.5 * delta_t,support,k2,params);  //calculate k2 at support and new t
	for(int i = 0; i < dimension; i++)
	{
		k2[i] *= delta_t; //multiply each element by delta_t
	}
	for(int i = 0; i < dimension; i++)
	{
		y[i] += k2[i]; 
	}
	free(support);
	free(k1);
	free(k2);
	



}





void rk4_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params)
{
	double a1 = 0.5;   //parameters
	double a2 = 0.5;
	double a3 = 1.0;
	double w1 = 1.0/6.0;
	double w2 = 1.0/3.0;
	double w3 = 1.0/3.0;
	double w4 = 1.0/6.0;
	double *support = malloc((int)(sizeof(double) * dimension));
	double *k1 = malloc((int)(sizeof(double) * dimension));
	double *k2 = malloc((int)(sizeof(double) * dimension));
	double *k3 = malloc((int)(sizeof(double) * dimension));
	double *k4 = malloc((int)(sizeof(double) * dimension));
	
	//#####################  k1 ###################
	
	func(t,y,k1,params);   //calculate k1
	for(int i = 0; i < dimension; i++)
	{
		k1[i] *= delta_t; //multiply each element by delta_t
	}
	
	
	//#####################  k2 ###################
	for(int i = 0; i < dimension; i++)
	
	{
		support[i] = y[i] + a1 * k1[i]; //fill support array for k2
	}
	func(t + a1 * delta_t,support,k2,params);  //calculate k2 at support and new t
	
	for(int i = 0; i < dimension; i++)
	{
		k2[i] *= delta_t; //multiply each element by delta_t
	}
	
	
	//#####################  k3 ###################
	
	
	for(int i = 0; i < dimension; i++)
	{
		support[i] = y[i] + a2 * k2[i]; //fill support array for k3
	}
	
	func(t + a2 * delta_t,support,k3,params);  //calculate k3 at support and new t
	
	for(int i = 0; i < dimension; i++)
	{
		k3[i] *= delta_t; //multiply each element by delta_t
	}
	
	
	//#####################  k4 ###################
	
	for(int i = 0; i < dimension; i++)
	{
		support[i] = y[i] + a3 * k3[i]; //fill support array for k4
	}
	
	func(t + a3 * delta_t,support,k4,params);  //calculate k4 at support and new t
	
	for(int i = 0; i < dimension; i++)
	{
		k4[i] *= delta_t; //multiply each element by delta_t
	}
	
	for(int i = 0; i < dimension; i++)
	{
		y[i] = y[i] +  w1 * k1[i] +  w2 * k2[i] + w3 * k3[i] + w4 * k4[i]; //calculate new y
	}
	free(support);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
	



}




void verlet_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params)
{
	
	int N = dimension/2; //no of x's and v's
	
	//acceleration arrays
	double *a1 = malloc((int)(sizeof(double) * dimension));
	double *a2 = malloc((int)(sizeof(double) * dimension));
	
	//fill a1
	func(t, y, a1, NULL);
	
	
	//calculate positions
	for(int i = 0; i < N; i++)
	{
		y[i] = y[i] + y[i+N] * delta_t + 0.5 * a1[i+N] * delta_t * delta_t;
	
	}
	
	//fill a2
	func(t, y, a2, NULL);
	
	//compute velocities
	for(int i = 0; i < N; i++)
	{
		y[i+N] = y[i+N] + 0.5  * delta_t * (a2[i+N] + a1[i+N]);
	
	}
	free(a1);
	free(a2);




}


int gauss(int n, int m, double A[n * m], int pivoting)
{	

	int sign = 1; //only for partial pivoting
	double c;
	double temp;
	//without partial pivoting
	if(pivoting == 0)
	{
		for(int i = 0; i < n - 1; i++)
		{
			for(int k = i + 1; k < n; k++)
			{
				c = A[k * m + i] / A[i * m + i];
				for(int j = 0; j < m; j++)
				{
					A[k * m + j] = A[k * m + j] - c * A[i * m + j];
				}
		
			}
	
	
		}
	return sign;
	}
	
	//with partial pivoting
	else
	{
		//rearrange matrix
		for(int i = 0; i < n - 1; i++)
		{
			for(int k = i + 1; k < n; k++)
			{
				if(fabs(A[i*m+i]) < fabs(A[k*m+i]))
				{	
					sign *= -1;
					for(int j = 0; j < m; j++)
					{
						temp = A[i * m + j];
						A[i * m + j] = A[k * m + j];
						A[k * m + j] = temp;
					}
		
				}
	
	
			}
	
		}
		
		//use gauss alg
		for(int i = 0; i < n - 1; i++)
		{
			for(int k = i + 1; k < n; k++)
			{
				c = A[k * m + i] / A[i * m + i];
				for(int j = 0; j < m; j++)
				{
					A[k * m + j] = A[k * m + j] - c * A[i * m + j];
				}
		
			}
	
	
		}


	return sign;
	}

}




double mat_det(int n, double A[n*n], int pivoting)
{
	int sign = 0;
	double det = 1;
	double m[n*n]; //placeholder matrix
	for(int i = 0; i < n * n; i++)
	{
		m[i] = A[i]; //copy matrix elements into m
	} 
	sign = gauss(n, n, m, pivoting);
	
	
	
	//compute determinant
	for(int i = 0; i < n; i++)
	{
		det *= m[i * n + i];
	}
	det *= sign;
	return det;







}



void lgs_solve(int n, double A[n*n], double b[n], double x[n])
{
	double C[n * (n+1)];
	double sum = 0.0;
	//first copy A into the nxn part of C
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			C[i * (n+1) + j] = A[i * n + j];
		}
	}
	
	//now add the b vector to the (n+1)th row
	for(int i = 0; i < n; i++)
	{
		C[i * (n + 1) + n] = b[i];
	}
	
	//triangularize
	gauss(n,n+1,C,1);
	
	
	
	
	//compute x elements
	x[n-1] = C[(n - 1) * (n+1) + n] / C[(n-1) * (n+1) + n - 1];
	for(int i = n-2; i >= 0; i--)
	{
		for(int j = i + 1; j < n; j++)
		{
			sum += C[i * (n+1) + j] * x[j];
		}
		x[i] = 1 / C[i * (n+1) + i] * (C[i * (n+1) + n] - sum);
		
		sum = 0;
	}
	

	return;



}


void mat_inverse(int n, double A[n*n], double A_inverse[n*n])
{
	double det = mat_det(n, A, 1);
	if(fabs(det) < 1e-5)
	{
		printf("Matrix not invertible!\n");
		return;
	}
	double b[n];
	double x[n];
	for(int i = 0; i < n ; i++)
	{
		b[i] = 0;
		x[i] = 0;
	}
	for(int i = 0; i < n ; i++)
	{
		b[i] = 1.0;
		lgs_solve(n, A, b, x);
		for(int j = 0; j < n; j++)
		{
			A_inverse[n*j + i] = x[j];
		}
		for(int k = 0; k < n ; k++)
		{
			b[k] = 0;
			x[k] = 0;
		}
	
	
	}
	return;
	




}




void lgs_jacobi_solve(int n, double A[], double b[], double x[], int k_max)
{
	int k = 0;
	double sum = 0; 
	double *xi = malloc((int)(sizeof(double) * n)); //array to store previous x in
	for(int i = 0; i < n; i++) xi[i] = x[i]; //copy x to xi
	while(k < k_max)
	{
	
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				if(i != j) sum += A[n * i + j] * xi[j]; //sum over i != j
			}
			x[i] = 1/A[i * n + i] * (b[i] - sum); //update x
			sum = 0;
		
		}
		for(int i = 0; i < n; i++) xi[i] = x[i]; //copy x to xi
		k++;
	}
	free(xi);
	return;




}



void lgs_gs_solve(int n, double A[], double b[], double x[], int k_max)
{
	int k = 0;
	double sum = 0; 
	while(k < k_max)
	{
	
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				if(i != j) sum += A[n * i + j] * x[j]; //sum over i != j
			}
			x[i] = 1/A[i * n + i] * (b[i] - sum); //update x
			sum = 0;
		
		}
		k++;
	}
	return;



}

