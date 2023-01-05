#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"
#include <math.h>





//constants
const double L = 5; //length
const int number = 6; //no of overall computation cycles
const double PI = 3.14159265358979323846;
const int dim = 500; //grid dimension
const double g = 0.5;



//Function to compute eigenvectors of H0 
void eigenvec(double A[], double x[], double w[]);


//Function to compute and print first four eigenvalues
void eigenval(double H0[], double x[], int N, FILE* ew);






//-------------------------------- MAIN--------------------------


int main()
{

	//output files
	FILE *ew = fopen("ew.dat","w");
	FILE *N_file = fopen("N.dat","w");

	
	int * N = (int*)malloc(number * sizeof(int)); //for each number of eigenvectors used we need different matrix dimensions
	
	
	//number of eigenvectors used
	for(int i = 0; i < number; i++)
	{
		N[i] = (int)pow(2.0, (double)(i+2));
	} 
	
	
	
	//compute eigenvalues and eigenvectors of H0 as in 28.1
	
	double dx = 2*L/(double)dim; //x step
	
	//x grid
	double x[dim]; 
	for(int i = 0; i < dim ; i++) x[i] = i * dx - L; 
	
	
	//H0 matrix
	double *A = (double *)malloc(dim*dim*sizeof(double));
	
	
	//eigenvalue array in ascending order
	double w[dim]; 
	
	//save Eigenvectors of H0 in A
	eigenvec(A,x,w);
	
	
	//compute eigenvalues of pertubed operator 
	
	
	for(int i = 0; i < number; i++)
	{
		eigenval(A,x,N[i],ew);
		fprintf(N_file, "%d\n",N[i]);
	}


	fclose(ew);
	fclose(N_file);
	free(N);
	free(A);
	


	return 0;
	
	
	
}





// ----------------------- EIGENVEC -----------------------------

void eigenvec(double A[], double x[], double w[])
{

	double dx = 2*L/(double)dim; //x step
	//lapack parameters
	lapack_int info,m,lda; 
	m = dim;
	lda = dim;
	
	
	// H0 in computational basis
	 
	for(int i = 0; i<m*m;i++) A[i] = 0; //set all values to zero
	for(int i = 0; i < m; i++) A[i*m+i] = 1 / (dx * dx) + 0.5 * x[i] * x[i]; //diagonal
	for(int i = 0; i < m-1; i++) A[i*m+i+1] = -0.5 / (dx*dx); //super diagonal
	
	
	//lapack subroutine
	
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'L', m, A, lda, w);
	return;

}



// ------------------------- 	EIGENVAL---------------------------------------------

void eigenval(double H0[], double x[], int N, FILE * ew)
{

	lapack_int info,m,lda; 
	m = N;
	lda = m;

	double mat_comp = 0; 
	double* H = (double*)malloc(sizeof(double) * N * N);
	double* w = (double*)malloc(sizeof(double) * N);
	for(int j = 0; j < N;j++)
	{
		for(int k = 0; k < N; k++)
		{
			//here we compute the jk matrix element
			for(int l = 0; l < dim; l++) mat_comp += H0[l*N + j] * g * x[l]* x[l]* x[l]* x[l] * H0[l*N + k];
			H[j*N+k] = mat_comp;
			//add diagonal term (eigenvalue of psi beta)
			if(j == k) H[j*N+k] += k + 0.5;
			mat_comp = 0;
		}
			
	}
	
	
	//lapack subroutine
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'U', m, H, lda, w);
		
	
	//print first 4 eigenvalues as a function of N
	
	for(int i = 0; i < 4; i++) fprintf(ew,"%lf\t",w[i]);
	fprintf(ew,"\n");
	free(H);
	free(w);
		
	
	return;



}

