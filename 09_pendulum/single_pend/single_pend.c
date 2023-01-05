#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <my_numerics.h>
#define PI 3.14159265358979323846


//constants

const int N = 1;              // no of pendulums
const double k = 100 ;           // spring constant
const double base_length = 0.0 ; // spring base length
const double mass = 1.0;        // pendulum mass
const double v0 = 1.0; //initial velocity



//output files
static const char euler_pos_file_name[] = "pos_euler.dat";
static const char rk2_pos_file_name[] = "pos_rk2.dat";
static const char rk4_pos_file_name[] = "pos_rk4.dat";





//physical system


int pendulums_ode(double t, const double y[], double f[], void *params)
{

    
 
	for(int i = 0; i < 2 * N; i++)
	{
		f[i] = 0; //initialize array
	}
    
	for(int i = 0; i < N; i++)
	{
		f[i] += y[N+i]; //fill velocities
	}
    
	f[N] += -k*y[0]/mass + k * (y[1] - y[0] - base_length)/ mass; //force on first mass 
    
	for(int i = 1; i < N - 1; i++)
	{
		f[N + i] +=  -k * (y[i] - y[i - 1] - base_length) / mass + k* (y[i+1] - y[i] - base_length) / mass; //fill forces
	}
	f[2 * N - 1] = -k * (y[N - 1]- y[N - 2] - base_length) / mass; //last pendulum is loose
	return 0;
}




//function to initialize arrays

void init(double y[], int N)
{
	
	for(int i = 0; i < N; i++)
	{
		y[i] = i * base_length;
	}
	for(int i = N + 1; i < 2 * N; i++)
	{
		y[i] = 0;

	}
    
	y[N] = v0;
	return;
}

// _---------------------- MAIN ------------------------------

int main(void)
{


	double T = 2 * PI * sqrt(mass / k); //period
	double Tmax = T/4.0; //fourth period
	double xmax = v0 * sqrt(mass / k); //x after T/4
	double delta_t[100];  //array for time steps
	delta_t[0] = 1;
	double step = pow(1.0e-8/delta_t[0],1./99.);
	for(int j = 0; j < 99; j++) delta_t[j + 1] = delta_t[j] * step;
	
	

 
	/*dimensionality*/
	int dimension = 2 * N;

  
    
	/* system arrays for all methods*/
	double y_eu[dimension];
	double y_rk2[dimension];
	double y_rk4[dimension];
    
	
    
	
  

  
	//open output files

	//open position output files
	FILE *pos_file_euler = fopen(euler_pos_file_name, "w");
	FILE *pos_file_rk2 = fopen(rk2_pos_file_name, "w");
	FILE *pos_file_rk4 = fopen(rk4_pos_file_name, "w");
    

	
	double t = 0.0; //simulation time
	double current_step = 0.0; //current delta_t
	double res_eu = 0.0; 
	double res_rk2 = 0.0; 
	double res_rk4 = 0.0; 
	int max_iter = 100; 
	
	
	
	for(int j = 0; j < max_iter; j++)
	{
		current_step = delta_t[j];
		
		//initialize
		init(y_eu,N);
		init(y_rk2,N);
		init(y_rk4,N);
	
		t = 0.0;
		
		int boool = 0;
		
		//actuals steps
		while (t <= Tmax)
		{
		
			if (boool) 
			{
                		res_eu = fabsf(y_eu[0] - xmax);
                		res_rk2 = fabsf(y_rk2[0] - xmax);
                		res_rk4 = fabsf(y_rk4[0] - xmax);
                		break;
           		}
           		
           		if ((Tmax - t) < delta_t[j])
           		{
               		 current_step = Tmax - t;
               		 boool = 1;
          		}	
			euler_step(t, current_step, y_eu, pendulums_ode, dimension, NULL);
			rk2_step(t, current_step, y_rk2, pendulums_ode, dimension, NULL);
			rk4_step(t, current_step, y_rk4, pendulums_ode, dimension, NULL);
			t += current_step;
		}
		
		fprintf(pos_file_euler, "%g\t%g\n", delta_t[j], res_eu);
		fprintf(pos_file_rk2, "%g\t%g\n", delta_t[j], res_rk2);
		fprintf(pos_file_rk4, "%g\t%g\n", delta_t[j], res_rk4);
		
		
	}
	

	fclose(pos_file_euler);
	fclose(pos_file_rk2);
	fclose(pos_file_rk4);
    	return 0;
  


}



