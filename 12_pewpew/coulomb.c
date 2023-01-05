#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <my_numerics.h>
#define PI 3.14159265358979323846


//constants

const double eps0 = 8.85419 * 1e-12; //vacuum permitivity
const double e = 1.602 * 1e-19; //elementary charge
const double m = 9.109e-31; //electron mass
const double delta_t = 1e-2;
const int N = 2; //dimension of ode system
double r0 = 10; //starting position
double rf = 10; //final position
double t0 = 0.0; //starting time
double tf = 2.0; //final time

//output files
static const char G_file_name[] = "G.dat";
static const char traj1_file_name[] = "traj1.dat";
static const char traj2_file_name[] = "traj2.dat";
static const char traj3_file_name[] = "traj3.dat";




//coulomb force

double F(double r)
{
	return e * e / (4 * PI * eps0 * r * r);
}

//physical system
int coulomb_ode(double t, const double y[], double f[], void *params)
{

	f[0] += y[1]; //compute velocity
	f[1] = 1 / m * F(y[0]); //compute force
	return 0;
}

//function to return difference y1(tf)-rf
double G(double v0)
{
	double t = t0; //system time
	double y[2] = {r0, v0}; //system array
	while(t < tf)
	{
		verlet_step(t, delta_t, y, coulomb_ode, N, NULL);
		t += delta_t;
	}
	return y[0] - rf;
}



int main()
{
	//set boundaries for v
	double vmin = -12.0; 
	double vmax = -2.0;
	double v = vmin;
	double delta_v = 1e-2;
	
	
	double x = 0.0; //will later be zero of G
	
	//compute and print zero of G
	x = find_root(G,-5.0,delta_t,1e-2,1000); 
	printf("%g \n",x);
	
	//open output files
	FILE *G_file = fopen(G_file_name, "w");
	FILE* traj[3] = {fopen(traj1_file_name, "w") , fopen(traj2_file_name, "w") , fopen(traj3_file_name, "w")};
	
	
	//compute G(v)
	while(v < vmax)
	{
		fprintf(G_file, "%g\t%g\n", v, G(v));
		v += delta_v;
	}
	
	double v0[3] = {-5.0,-4.0,x};
	//define system array to compute trajectories
	double y[2];
	double t = 0.0;
	
	
	for(int i = 0; i < 3; i++)
	{
		//intialize y
		y[0] = r0;
		y[1] = v0[i];
		
		//compute trajectories
		while(t <= tf)
		{	
			fprintf(traj[i], "%g\t%g\n", t, y[0]);
			verlet_step(t, delta_t, y, coulomb_ode, N, NULL);
			t += delta_t;
			
		}
	fclose(traj[i]);
	t = 0.0;
	
	}
	
	
	
	return 0;
	


}
