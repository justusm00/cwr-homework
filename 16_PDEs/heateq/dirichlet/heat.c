#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>

//constants 
const double L = 1; //rod length
const double D = 0.1; //heat constant
const double dx = 0.01; //x step
const int N = (int)L/dx; //no of points
const double Tl = 1; //const temperature at left boundary
const double Tr = -1; //const temperature at right boundary
const double tmax = 1; //maximum simulation time


void heat_ftcs(double t, double y[], double dt)
{
	double temp[N];
	temp[0] = Tl; //first compute left value 
	for(int i = 1; i < N - 1; i++)
	{
		temp[i] = y[i] + dt * D * (y[i - 1] - 2 * y[i] + y[i+1]) / (dx * dx); //compute every other value except the rightmost
	
	}
	temp[N-1] = Tr; //rightmost 
	for(int i = 0; i < N; i++) y[i] = temp[i]; //overwrite y
	return;

}

int main()
{
	double y[N];
	int j = 0; //iteration counter
	//initialize array
	for(int i = 1; i < N - 1; i++) y[i] = 0; 
	y[0] = Tl;
	y[N-1] = Tr;
	
	
	double dt[2]; //array for two different time steps
	dt[0] = 0.8 * dx * dx / (2 * D); //stable
	dt[1] = 1.01 * dx * dx / (2 * D); //unstable
	printf("time steps: %lf\t%lf\n",dt[0],dt[1]);
	double t = 0; //system time
	FILE* heatfile[2] = {fopen("heat1.dat", "w") , fopen("heat2.dat", "w")}; //output files
	FILE* timefile[2] = {fopen("time1.dat", "w") , fopen("time2.dat", "w")}; //output files
	
	//first do stable simulation
	while(t<tmax)
	{
		if((j % 10) == 0) //only print a few values as otherwise we get way too many data points
		{
			for(int i = 0; i < N-1; i++) fprintf(heatfile[0],"%lf,",y[i]);
			fprintf(heatfile[0],"%lf",y[N-1]);
			fprintf(heatfile[0],"\n");
			fprintf(timefile[0],"%lf\n",t);
		}
		heat_ftcs(t,y,dt[0]);
		t += dt[0];
		j++;
	}
	
	
	
	//reinitialize everything and repeat for unstable dt
	
	//initialize array
	for(int i = 1; i < N - 1; i++) y[i] = 0; 
	y[0] = Tl;
	y[N-1] = Tr;
	t = 0;
	j = 0;
	
	while(t<tmax)
	{
		if((j % 10) == 0) //only print a few values as otherwise we get way too many data points
		{
			for(int i = 0; i < N-1; i++) fprintf(heatfile[1],"%lf,",y[i]);
			fprintf(heatfile[1],"%lf",y[N-1]);
			fprintf(heatfile[1],"\n");
			fprintf(timefile[1],"%lf\n",t);
		}
		
		heat_ftcs(t,y,dt[1]);
		t += dt[1];
		j++;
	}
	fclose(heatfile[0]);
	fclose(heatfile[1]);
	fclose(timefile[0]);
	fclose(timefile[1]);
	return 0;


}
