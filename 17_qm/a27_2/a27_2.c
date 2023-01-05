#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"


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

double w_anal(n)
{
	return 0.5 * PI * PI * n * n;
}

int main()
{
	FILE *ew = fopen("ew.dat","w");
	FILE *ev = fopen("ev.dat","w");
	int N = 1000; //dimension
	lapack_int info,m=N,lda=m; //lapack parameters
	double w[N]; //eigenvalues in ascending order#
	
	double dx = L/((double)N); //x step
	
	double *A = (double *)malloc(m*m*sizeof(double)); //matrix
	
	
	
	//fill matrix
	for(int i = 0; i<N*N;i++) A[i] = 0; //set all values to zero
	for(int i = 0; i < N; i++) A[i*N+i] = 1 / (dx * dx); //diagonal
	for(int i = 0; i < N-1; i++) A[i*N+i+1] = -0.5 / (dx*dx); //super diagonal
	
	/*
	print_matrix_rowmajor("Entry Matrix A:",m,m,A,m);
	*/
	
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', m, A, lda, w);
	
	

	
	//print first 20 eigenvalues and the according analytical solution
	for(int i = 0; i < 20; i++) fprintf(ew,"%lf\t%lf\n",w[i],w_anal(i+1));
	
	//print first 3 eigenvectors
	for(int i = 0; i < N; i++)
	{
		fprintf(ev,"%lf\t%lf\t%lf\n",A[i*N], A[i*N + 1], A[i*N + 2]);
	}
	fclose(ew);
	fclose(ev);
	return 0;


}
