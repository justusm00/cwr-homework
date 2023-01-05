#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <lapacke.h>



/*-----------------------------------------------------------------------------*/


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


double solve_linear_system(lapack_int n)
{

	clock_t begin, end; //time variables
	double cpu_time;
	begin = clock();
	// Variables needed for LAPACK
	// - info: takes the returned value from the LAPACK subroutine
	// - m: number of linear equations  
	// - lda: The leading dimension of the matrix A
	// - ldb: The leading dimension of the array b
	// - nrhs: the number of columns of the right hand side
	lapack_int info,m=n,lda=m,ldb=1,nrhs=1;
	double *A = (double *)malloc(m*m*sizeof(double));
 	double *b = (double *)malloc(m*sizeof(double));
	FILE *values = fopen("matrices.dat","r");
	
	
	// ipiv is an array of dimension m that stores what lines were
	// changed during the pivoting procedure
	lapack_int *ipiv = (lapack_int *)malloc(m*sizeof(lapack_int)) ;
	
	for(int i = 0; i < m; i++)
	{
		for(int j = 0; j < m; j++) fscanf(values, "%lf", &A[i*m+j]); //fill arrays with random values
		fscanf(values, "%lf", &b[i]);
	}
	

	/*
	print_matrix_rowmajor("Entry Matrix A:",m,m,A,m);
	print_matrix_rowmajor("vector b:",m,1,b,1);
	*/
	
	// The majority of the arguments here where explained at the beginning 
	// of main(), the only left are the first the second parameters, the former 
	// specify how the input vectors (A and b) are organized. Row major means
	// that we have inserted values in A and b going through the rows, i.e. the
	// first m elements of A are the firts row of the matrix. The second arument
	// 'U' says to LAPACK to store the matrix A in the upper form. In this case 
	// it doesn't change the result.
	info =  LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', m, nrhs, A, lda, ipiv, b, ldb);
  	/*
	print_matrix_rowmajor( "Solution: ", m, 1, b, 1 );
	*/
  
	fclose(values);
	free(A);
	free(b);
	free(ipiv);
 
	
	
	
	
	
	end = clock();
	cpu_time = (double)(end - begin)/CLOCKS_PER_SEC; //computation time
	return cpu_time;



}


int main()
{


	FILE *time = fopen("cpu_time.dat","w");
	double t = 0;
	int n = 500; //dimension
	int number = 10; //no of calculations
	int step = 450;
	t = solve_linear_system(n);
	
	
	for(int i = 0; i < number; i++)
	{
		
		t = solve_linear_system(n);
		
		printf("%lf\n",t);
		fprintf(time,"%d\t%lf\n",n,t);
		n += step;
		
	}
	
	fclose(time);
	return 0;


}
