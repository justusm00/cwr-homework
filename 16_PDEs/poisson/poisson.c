#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>

//constants
const double l_el = 0.002; //electrode dimension in cm
const double l_box = 0.01; //box dimension in cm
const double phi_el = 1; //electrode potential in volts
const double phi_box = 0; //box potential in volts

const double eps0 = 8.854e-12; //vacuum permeability




//function to do gauss seidel iteration
void poisson_gs(double dx, double dy, double phi[], int n, double el_min, double el_max)
{

	double c = 0.5 * 1 / (dx * dx + dy * dy); //prefactor
	for(int i = 1; i < n-1; i++)
	{
		for(int j = 1; j < n-1; j++)
		{
			if((i > el_min-1) && (j > el_min-1) && (i < el_max-1) && (j < el_max-1))
			{}
			else phi[i*n + j] = c * (dx * dx * (phi[(i+1)*n+j] + phi[(i-1)*n+j]) + dy * dy * (phi[i*n+j+1] + phi[i*n+j-1]));
		
		}
	}
	return;


}


int main()

{
	double dx = 1e-4; //x step
	double dy = 1e-4; //y step
	double max_iter = 5000;
	double rel_tol = 1e-5; //stop iteration if residual becomes smaller
	int k = 0; //iteration counter
	int n = (int)(l_box/dx); //grid dimension
	int n_el = (int)(l_el/dx); //electrode dimension
	double *phi = (double*)malloc(n*n*sizeof(double)); //potential array
	double *phi_old = (double*)malloc(n*n*sizeof(double)); //previous step to compute modulus
	double *ex = (double*)malloc(n*n*sizeof(double)); //x component of electric field
	double *ey = (double*)malloc(n*n*sizeof(double)); //y component of electric field
	double modulus = 0;
	FILE* pot_file = fopen("pot.dat","w");
	FILE* ex_file = fopen("ex.dat","w");
	FILE* ey_file = fopen("ey.dat","w");
	FILE* e_file = fopen("e.dat","w"); //abs value of E
	
	//define index limits to iterate through the electrode
	int el_min = (int)(n - n_el)/2; //minimal electrode index
	int el_max = n - el_min; //maximal electrode index
	printf("%d\n%d\n",el_min,el_max);
	
	//first set all values of phi to zero
	//this way, the box boundaries are already correct
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
		
			if((i > el_min-1) && (j > el_min-1) && (i < el_max-1) && (j < el_max-1))
			{
			
				phi[i*n + j] = 1;
			}
			else phi[i*n + j] = 0;
		}
	}
	
	
	//do first step outside the loop to compute first modulus
	for(int i = 0; i < n*n; i++)
	{
		phi_old[i] = phi[i];
	}
		
	poisson_gs(dx,dy,phi,n,el_min,el_max);
	for(int i = 0; i < n*n; i++)
	{
		modulus += (phi_old[i]-phi[i])*(phi_old[i]-phi[i]); //compute modulus
	}
	modulus = sqrt(modulus);
	
	
	
	//do iterations 
	while((k < max_iter) && (modulus >= rel_tol))
	{
		modulus = 0;
		for(int i = 0; i < n*n; i++)
		{
			phi_old[i] = phi[i];
		}
		
		poisson_gs(dx,dy,phi,n,el_min,el_max);
		for(int i = 0; i < n*n; i++)
		{
			modulus += (phi_old[i]-phi[i])*(phi_old[i]-phi[i]); //compute modulus
		}
		modulus = sqrt(modulus);
		k++;
	
	}
	
	printf("No of iterations: %d\n",k);
	printf("final modulus: %lf\n",modulus);

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n-1; j++)
		{
			fprintf(pot_file,"%lf,",phi[i*n + j]);
		}
		fprintf(pot_file,"%lf",phi[i*n + n - 1]);
		fprintf(pot_file,"\n");
	}
	
	//initialize field arrays
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			ex[i * n + j] = 0;
			ey[i * n + j] = 0;
		}
	}
	
	//compute electric field
	for(int i = 1; i < n-1; i++)
	{
		for(int j = 1; j < n-1; j++)
		{
			ey[i * n + j] = -(phi[(i+1)*n + j] - phi[(i-1)*n + j]) / (2*dx);
			ex[i * n + j] = -(phi[i*n + j + 1] - phi[i*n + j - 1]) / (2*dy);
		}
	}
	
	
	//print field
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n-1; j++)
		{
			fprintf(ex_file,"%lf,",ex[i*n + j]);
			fprintf(ey_file,"%lf,",ey[i*n + j]);
			fprintf(e_file,"%lf,",sqrt(ex[i*n + j] * ex[i*n + j] + ey[i*n + j] * ey[i*n + j]));
		}
		fprintf(ex_file,"%lf",ex[i*n + n - 1]);
		fprintf(ey_file,"%lf",ey[i*n + n - 1]);
		fprintf(e_file,"%lf", sqrt(ex[i*n + n - 1] * ex[i*n + n - 1] + ey[i*n + n - 1] * ey[i*n + n - 1]));
		fprintf(ex_file,"\n");
		fprintf(e_file,"\n");
		fprintf(ey_file,"\n");
	}
	
	fclose(pot_file);
	fclose(ex_file);
	fclose(ey_file);
	fclose(e_file);
	free(phi);
	free(phi_old);
	free(ex);
	free(ey);
	return 0;
}
