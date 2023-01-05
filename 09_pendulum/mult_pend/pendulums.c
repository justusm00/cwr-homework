#include <stdio.h>
#include <stdlib.h>
#include <my_numerics.h>


//constants

const int N = 30;              // no of pendulums
const double k = 100 ;           // spring constant
const double base_length = 1 ; // spring base length
const double mass = 1;        // pendulum mass

//parameters

const double T_max = 20.0;
const double delta_t = 1e-2; //time step


//output files
static const char euler_pos_file_name[] = "euler_pos.dat";
static const char rk2_pos_file_name[] = "rk2_pos.dat";
static const char rk4_pos_file_name[] = "rk4_pos.dat";
static const char euler_energy_file_name[] = "euler_energy.dat";
static const char rk2_energy_file_name[] = "rk2_energy.dat";
static const char rk4_energy_file_name[] = "rk4_energy.dat";
static const char time_file_name[] = "time.dat";

static const char step_file_name[] = "delta.dat"; //print step just for fun



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


//energy function


double pendulums_energy(const double y[])
{
	double energy = 0;
	double v = 0;

	energy += 0.5 * k * y[0] * y[0];
	for(int i = 0; i < N; i++)
      	  {
        		v = y[i+N];
        		energy += 0.5 * mass * v * v; //kinetic energy
        		energy += 0.5 * k * (y[i + 1]- y[i] - base_length) * (y[i + 1]- y[i] - base_length); //potential energy
      	 }
	return energy;
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
    
	y[N] = 20;
	return;
}

// _---------------------- MAIN ------------------------------

int main(void)
{
    

 
	/*dimensionality*/
	int dimension = 2 * N;

  

  

    
	/* system arrays for all methods*/
	double y_eu[dimension];
	double y_rk2[dimension];
	double y_rk4[dimension];
    
	//initial values
    
	init(y_eu,N);
	init(y_rk2,N);
	init(y_rk4,N);
  

   
    
    
	//open output files

	//open position output files
	FILE *pos_file_euler = fopen(euler_pos_file_name, "w");
	FILE *pos_file_rk2 = fopen(rk2_pos_file_name, "w");
	FILE *pos_file_rk4 = fopen(rk4_pos_file_name, "w");
    
	//open energy output files
	FILE *energy_file_euler = fopen(euler_energy_file_name, "w");
	FILE *energy_file_rk2 = fopen(rk2_energy_file_name, "w");
	FILE *energy_file_rk4 = fopen(rk4_energy_file_name, "w");
    
	//open time file
	FILE *time_file = fopen(time_file_name, "w");
	FILE *step_file = fopen(step_file_name, "w");
    
	fprintf(step_file, "%g", delta_t);
	

	// simulation time
	double t = 0;
  

	//actual steps
	while (t < T_max)
	{
		//do steps
		euler_step(t, delta_t, y_eu, pendulums_ode, dimension, NULL);
		rk2_step(t, delta_t, y_rk2, pendulums_ode, dimension, NULL);
		rk4_step(t, delta_t, y_rk4, pendulums_ode, dimension, NULL);

	     
		//print position
		for(int i = 0; i < N; i++)
        	{

			fprintf(pos_file_euler, "%g \t", y_eu[i]);
        		fprintf(pos_file_rk2, "%g \t",y_rk2[i]);
        		fprintf(pos_file_rk4, "%g  \t",y_rk4[i]);
        	}
        
        	//print energy
        	fprintf(energy_file_euler, "%g \t %g\n", t, pendulums_energy(y_eu));
        	fprintf(energy_file_rk2, "%g \t %g\n", t, pendulums_energy(y_rk2));
        	fprintf(energy_file_rk4, "%g \t %g\n", t, pendulums_energy(y_rk4));
        
        	//make new line
        	fprintf(pos_file_euler, "\n");
        	fprintf(pos_file_rk2, "\n");
        	fprintf(pos_file_rk4, "\n");
        
        	//print time in extra file
        	fprintf(time_file, "%g \n",t);
        
        	t += delta_t;
	}

	fclose(pos_file_euler);
	fclose(pos_file_rk2);
	fclose(pos_file_rk4);
	fclose(energy_file_euler);
	fclose(energy_file_rk2);
	fclose(energy_file_rk4);
	fclose(time_file);
	fclose(step_file);
	return 0;
}



