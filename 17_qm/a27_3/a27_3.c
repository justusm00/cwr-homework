#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"
#include <math.h>


//constants
double L = 1; //length
const double PI = 3.14159265358979323846;

// Just a simple function to print out a matrix/vector

void print_matrix_rowmajor(char *desc, lapack_int m, lapack_int n, double *mat,
                           lapack_int ldm)	
{
	lapack_int i, j;
	printf( "\n %s\n", desc ); 
	for( i = 0; i < m; i++ ) 
	{
		for( j = 0; j < n; j++) printf( " %6.2f", mat[i*ldm+j]);
		printf( "\n" );
	}
}

double w_anal(int n)
{
	return 0.5 * PI * PI * n * n;
}

int main()
{
	int number = 10; //number of dxs
	int N[number]; //dimension array to compute eigenvalues for ten different dx
	double dx = 1e-3; //x step, will later be varied
	for(int i = 0; i < number; i++)
	{
		N[i] = (int)(L/dx);
		dx += 0.009;
	}
	
	//create individual eigenvalue-arrays for each N
	double w0[N[0]];
	double w1[N[1]];
	double w2[N[2]];
	double w3[N[3]];
	double w4[N[4]];
	double w5[N[5]];
	double w6[N[6]];
	double w7[N[7]];
	double w8[N[8]];
	double w9[N[9]];

	double *w[] = {w0, w1,w2,w3,w4,w5,w6,w7,w8,w9};
	
	FILE *ew = fopen("ew.dat","w");
	FILE *dx_file = fopen("dx.dat","w");
	lapack_int info,m,lda; //lapack parameters
	
	
	//create individual matrices for each N 
	double *A0 = (double *)malloc(N[0]*N[0]*sizeof(double)); 
	double *A1 = (double *)malloc(N[1]*N[1]*sizeof(double)); 
	double *A2 = (double *)malloc(N[2]*N[2]*sizeof(double)); 
	double *A3 = (double *)malloc(N[3]*N[3]*sizeof(double)); 
	double *A4 = (double *)malloc(N[4]*N[4]*sizeof(double)); 
	double *A5 = (double *)malloc(N[5]*N[5]*sizeof(double)); 
	double *A6 = (double *)malloc(N[6]*N[6]*sizeof(double)); 
	double *A7 = (double *)malloc(N[7]*N[7]*sizeof(double)); 
	double *A8 = (double *)malloc(N[8]*N[8]*sizeof(double)); 
	double *A9 = (double *)malloc(N[9]*N[9]*sizeof(double)); 
	
	double *A[] = {A0,A1,A2,A3,A4,A5,A6,A7,A8,A9};
	
	
	//main loop
	for(int k = 0; k < number; k++)
	{
		m = N[k];
		lda = m;
		dx = L/((double)m);
		//fill matrix
		for(int i = 0; i<m*m;i++) A[k][i] = 0; //set all values to zero
		for(int i = 0; i < m; i++) A[k][i*m+i] = 1 / (dx * dx); //diagonal
		for(int i = 0; i < m-1; i++) A[k][i*m+i+1] = -0.5 / (dx*dx); //super diagonal
	

		info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'L', m, A[k], lda, w[k]);
	
	
		//print first 4 errors  (numerical eigenvalue/ analytical eigenvalue)
		fprintf(dx_file,"%lf\n",dx);
		for(int i = 0; i < 4; i++) fprintf(ew,"%lf\t", w[k][i]/(w_anal(i+1)));
		fprintf(ew,"\n");
	
	}
	fclose(ew);
	fclose(dx_file);
	return 0;


}
