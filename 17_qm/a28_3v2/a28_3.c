#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "lapacke.h"
#include <math.h>







//constants
const double L = 5.0; //length
const int number = 6; //no of overall computation cycles
const double PI = 3.14159265358979323846;
const int dim = 500; //grid dimension
const double g = 0.5;



//Function to compute eigenvectors of H0 
void eigenvec(double A[], double x[], double w[]);


//Function to fill H (operator with disturbance) 
void eigenval(double H0[], double x[], int N, FILE* ew);






//-------------------------------- MAIN--------------------------


int main()
{

	//output files
	FILE *ew = fopen("ew.dat","w");
	FILE *ev = fopen("ev.dat","w");
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
	
	
	eigenvec(A,x,w);
	
	
	
	//---------------------use eigevectors of H0 to get H -------------------------------------------------
	
	
	
	for(int k = 0; k < number; k++)
	{
		
		eigenval(A, x, N[k], ew);
		fprintf(N_file, "%d\n",N[k]);
			
	}
	
	
	fclose(ew);
	fclose(ev);
	fclose(N_file);
	
	
	free(A);
	free(N);
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



// ------------------------- EV---------------------------------------------

void eigenval(double H0[], double x[], int N, FILE* ew)
{

	double dx = 2*L/(double)N; //x step
	double mat_comp = 0;
	double bla = 0;
	lapack_int info,m,lda; 
	m = N;
	lda = N;
	double * Hu = (double*)malloc(sizeof(double) * dim);; //compontents of H|ui>
	double * H = (double*)malloc(sizeof(double) * dim * dim); //array for H
	double * H1 = (double*)malloc(sizeof(double) * N * N); //array for H in eigenbasis of H0
	double * w = (double*)malloc(sizeof(double) * N); //array for w
	

	//fill H 
	for(int i = 0; i<dim*dim;i++) H[i] = 0; //set all values to zero
	for(int i = 0; i < dim; i++) H[i*dim+i] = 1 / (dx * dx) + 0.5 * x[i] * x[i] + g * x[i] * x[i] * x[i] * x[i]; //diagonal
	for(int i = 0; i < dim-1; i++) H[i*dim+i+1] = -0.5 / (dx*dx); //super diagonal
	for(int i = 1; i < dim; i++) H[i*dim+i-1] = -0.5 / (dx*dx); //sub diagonal
	
	
	
	//get H into basis of H0
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
		
			//compute H|u_j>
			for(int k = 0; k < dim; k++)
			{
				for(int l = 0; l < dim; l++)
				{
					bla += H[k*dim+l] * H0[l*dim+j];
				}
				Hu[k] = bla;
				bla = 0;
			}
			//here we compute the ij matrix element
			for(int k = 0; k < dim; k++) mat_comp += H0[k*dim+i] * Hu[k];
			H1[i*N+j] = mat_comp;
			mat_comp = 0;
		}
	}
	
	
	
	
	
	//compute eigenvalues with lapack
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'N', 'L', m, H1, lda, w);
	//print first 4 eigenvalues as a function of N
	
	for(int i = 0; i < 4; i++) fprintf(ew,"%lf\t",w[i]);
	
	
	fprintf(ew,"\n");
	free(H);
	free(H1);
	free(w);
	free(Hu);
	return;



}
