#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <my_numerics.h>



//constants
const int N = 4; //system dimension
const double PI = 3.14159265358979323846;
const double g = 9.81; //gravitational acceleration
const double cw = 0.47; //resistance value for even sphere
const double rho = 1.2041; //density of air
const double R = 0.055; //cannon ball radius
const double v0 = 250.0; //initial velocity
const double m = 5.5; //cannon ball mass
const double Re = 1e6; //reynolds number


//output files
static const char angle_file_name[] = "angles.dat";
static const char pos_file_name[] = "pos.dat";




//define absolute value of 2D vector
double vecabs(double x,double y)
{
	double absval = sqrt(x * x + y * y);
	return absval;
}

//system function

int cannonball_ode(double t, const double y[], double f[], void *params)
{
	double v = vecabs(y[2],y[3]); //current velocity
	double A = PI * R * R; //area
	double Fd = rho * v * v / 2 * A * cw; //absolute value of friction force
	
	f[0] = y[2];
	f[1] = y[3]; //fill velocities
	f[2] = -Fd * y[2]/v; //just friction force in x direction
	f[3] = -m * g - Fd * y[3]/v; //gravitational force and friction force in y direction
	return 0;

}


int main()
{
	int boool = 1;
	double alpha[10]; //array of initial angles
	double vx = 0; //initial x-velocity
	double vy = 0; //initial y-velocity
	double y[10][N]; //10 trajectories
	double delta_t = 1e-3; //step
	double t = 0.0;
	
	
	
	//open output files
	FILE *pos_file = fopen(pos_file_name,"w");
	FILE *angle_file = fopen(angle_file_name,"w");
	
	//fill angles and print them
	for(int i = 0; i < 10; i++)
	{
		alpha[i] = (9 * i + 1)*PI/180 ; //angle in rad
		fprintf(angle_file,"%g\n", alpha[i]);
	}
	
	//initialize y's
	for(int i = 0; i < 10; i++)
	{
		vx = v0 * cos(alpha[i]);
		vy = v0 * sin(alpha[i]);
		y[i][0] = 0.0;
		y[i][1] = 0.0;
		y[i][2] = vx;
		y[i][3] = vy;
	}
	
	//do first step
	for(int i = 0; i < 10; i++)
	{
		rk2_step(t, delta_t, y[i], cannonball_ode, N, NULL);
	}
	
	
	
	
	
	
	while(boool < 10)
	{	
		boool = 0;
		//do step
		for(int i = 0; i < 10; i++)
		{
			if(y[i][1] > 1e-4) //if the current y value isnt zero, go one step ahead
			{
				rk2_step(t, delta_t, y[i], cannonball_ode, N, NULL);
			}	
			else boool++; //count one up, if all values are zero then the loop stops
			t += delta_t;
			fprintf(pos_file, "%g\t%g\t",y[i][0],y[i][1]);
			
		
		}
		
		fprintf(pos_file, "\n");
	

	
	}
	fclose(angle_file);
	fclose(pos_file);
	return 0;



}
