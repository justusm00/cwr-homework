#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"


//---------------- here I do the convergence stuff

struct arbitrary_len_arr {
    int len;
    double* arr;
};


//constants
const double L = 5.0; //length
const double PI = 3.14159265358979323846;
const double g = 0.5;





int main()
{
	int number = 1000; //number of dx's
	int N[number]; //dimension array to compute eigenvalues for ten different dx
	double dx = 1e-2; //x step, will later be varied
	double dxmax = 1;
	double dxs[number]; //array with different x steps
	double dx_step = (dxmax-dx)/number; //increment for dx
	for(int i = 0; i < number; i++)
	{
		dxs[i] = dx;
		N[i] = (int)2 * L / dx; //dimension that belongs to dx_i
		dx += dx_step;
	}
	
	
	
	//super big brain time :D
	struct arbitrary_len_arr* w = malloc(sizeof(struct arbitrary_len_arr) * number); //pointer array to arrays that the eigenvalues will be stored in
	struct arbitrary_len_arr* x = malloc(sizeof(struct arbitrary_len_arr) * number); //pointer array to arrays that the x grids will be stored in
	struct arbitrary_len_arr* A = malloc(sizeof(struct arbitrary_len_arr) * number); //pointer array to arrays that the matrices will be stored in
	
	
	//allocate memory for each array
	for(int i = 0; i < number; i++)
	{
		w[i].len = N[i];
		w[i].arr = malloc(sizeof(double) * N[i]);
		x[i].len = N[i];
		x[i].arr = malloc(sizeof(double) * N[i]);
		A[i].len = N[i] * N[i];
		A[i].arr = malloc(sizeof(double) * N[i] * N[i]);
	}
	
	
	
	//output files
	FILE *ew = fopen("ew.dat","w");
	FILE *dx_file = fopen("dx.dat","w");
	
	
	//lapack parameters
	lapack_int info,m,lda; 
	
	
	
	//main loop
	for(int k = 0; k < number; k++)
	{
		m = N[k];
		lda = m;
		dx = dxs[k];
		for(int i = 0; i < N[k] ; i++) x[k].arr[i] = i * dx - L;
		//fill matrix
		for(int i = 0; i<m*m;i++) A[k].arr[i] = 0; //set all values to zero
		for(int i = 0; i < m; i++) A[k].arr[i*m+i] = 1 / (dx * dx) + g * x[k].arr[i]* x[k].arr[i]* x[k].arr[i]* x[k].arr[i]; //diagonal
		for(int i = 0; i < m-1; i++) A[k].arr[i*m+i+1] = -0.5 / (dx*dx); //super diagonal
	

		info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'L', m, A[k].arr, lda, w[k].arr);
	
	
		//print the relative residues for neigbouring delta xi's 
		if(k > 0)
		{
			fprintf(dx_file,"%lf\n",dx);
			for(int i = 0; i < 4; i++) fprintf(ew,"%lf\t", fabsf((w[k].arr[i]-w[k-1].arr[i])/w[k].arr[i]));
			fprintf(ew,"\n");
		}
	
	}
	fclose(ew);
	fclose(dx_file);
	
	//free the arrays for they have served us well
	for(int i = 0; i < number; i++)
	{
		free(A[i].arr);
		free(x[i].arr);
		free(w[i].arr);
		
	}
	return 0;


}
