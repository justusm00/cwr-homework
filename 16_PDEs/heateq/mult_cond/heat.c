#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>

//constants 
const double L = 1; //rod length
const double D = 0.1; //heat constant
const double dx = 0.01; //x step
const int N = (int)L/dx; //no of points
const double tmax = 10; //maximum simulation time




double peak(double t)
{
	return exp(-(t - 5) * (t - 5));
}


void heat_ftcs(double t, double y[], double dt)
{
	double temp[N];
	temp[0] = peak(t); //dirichlet
	for(int i = 1; i < N - 1; i++)
	{
		temp[i] = y[i] + dt * D * (y[i - 1] - 2 * y[i] + y[i+1]) / (dx * dx); //compute every other value except the rightmost
	
	}
	temp[N-1] = temp[N-2]; //neumann 
	for(int i = 0; i < N; i++) y[i] = temp[i]; //overwrite y
	return;

}



int main()
{
	double y[N];
	int j = 0; //iteration counter
	//initialize array
	for(int i = 0; i < N ; i++) y[i] = 0; 
	
	

	
	double dt; //time step
	dt = 0.8 * dx * dx / (2 * D); //stable

	double t = 0; //system time
	FILE* heatfile = fopen("heat.dat", "w");
	FILE* timefile= fopen("time.dat", "w");
	
	//simulation
	while(t<tmax)
	{
		if((j % 100) == 0) //only print a few values as otherwise we get way too many data points
		{
			for(int i = 0; i < N-1; i++) fprintf(heatfile,"%lf,",y[i]);
			fprintf(heatfile,"%lf",y[N-1]);
			fprintf(heatfile,"\n");
			fprintf(timefile,"%lf\n",t);
		}
		heat_ftcs(t,y,dt);
		t += dt;
		j++;
	}
	
	
	
	
	fclose(heatfile);
	fclose(timefile);
	return 0;


}
