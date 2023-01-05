#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"


// ------------ exactly the same as 28.1 but with g = 0.5 -------------------

//constants
double L = 5; //length
const double PI = 3.14159265358979323846;
const int N = 500; //dimension
const double g = 0.5;




int main()

{
	double dx = 2*L/(double)N; //x step
	
	//x grid
	double x[N]; 
	for(int i = 0; i < N ; i++) x[i] = i * dx - L; 
	
	
	//lapack parameters
	lapack_int info,m=N,lda=m; 
	
	//hamilton operator in computational basis
	double *A = (double *)malloc(m*m*sizeof(double)); 
	for(int i = 0; i<m*m;i++) A[i] = 0; //set all values to zero
	for(int i = 0; i < m; i++) A[i*m+i] = 1 / (dx * dx) + 0.5 * x[i] * x[i] + g * x[i] * x[i] * x[i] * x[i]; //diagonal
	for(int i = 0; i < m-1; i++) A[i*m+i+1] = -0.5 / (dx*dx); //super diagonal
	
	double w[N]; //eigenvalues in ascending order
	
	//output files
	FILE *ew = fopen("ew.dat","w");
	FILE *ev = fopen("ev.dat","w");
	
	
	//lapack subroutine
	
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', m, A, lda, w);
	
	
	//print first 20 eigenvalues and the according analytical solution
	for(int i = 0; i < 20; i++) fprintf(ew,"%lf\n",w[i]);
	
	//print absolute value squared of first 3 eigenvectors
	for(int i = 0; i < N; i++)
	{
		fprintf(ev,"%lf\t%lf\t%lf\n",A[i*N] * A[i*N], A[i*N + 1] * A[i*N + 1], A[i*N + 2] * A[i*N + 2]);
	}
	fclose(ew);
	fclose(ev);
	return 0;

}
