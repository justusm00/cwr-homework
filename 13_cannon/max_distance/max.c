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
const double delta_t = 1e-1; //time step
const double delta = 1e-3; //diff step





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


//distance function

double cannon_distance(double angle)
{
	double xf = 0.0; //distance
	double vx = v0 * cos(angle); //initial x-velocity
	double vy = v0 * sin(angle);  //initial y-velocity
	double y[N];
	
	//initialize y
	y[0] = 0.0;
	y[1] = 0.0;
	y[2] = vx;
	y[3] = vy;
	double t = 0.0;
	//place holders for interpolation
	double xi;
	double yi;
	rk2_step(t, delta_t, y, cannonball_ode, N, NULL);
	while(y[1] > 1e-4)
	{
		rk2_step(t, delta_t, y, cannonball_ode, N, NULL);
		t += delta_t;
	}
	xi = y[0];
	yi = y[1];
	
	//one more step
	rk2_step(t, delta_t, y, cannonball_ode, N, NULL);
	//interpolate
	xf = (xi * y[1] - yi * y[0]) / (y[1] - yi);
	return xf;
}


//function to return derivative of the distance function
double D_cannon_distance(double angle)
{
	double y = 0.0;
	y = diff(angle, delta, cannon_distance);
	return y;

}

int main()
{
	double alpha_max = 0.0; //maximum angle
	double d_max = 0.0; //maximum distance
	double x0 = 35 * PI / 180;
	alpha_max =  find_root(D_cannon_distance,x0, delta, 1e-1,1000);
	d_max = cannon_distance(alpha_max);
	printf("Der maximale Winkel liegt bei %g ° , die maximale Strecke bei %g m \n", alpha_max * 180 / PI, d_max);
	return 0;


}
