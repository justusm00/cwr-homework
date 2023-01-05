#include<stdio.h>
#include<stdlib.h>
#include<tgmath.h>


const double a = 1.0/6.1;
const double gamma = 1.0/3.7;
const double N = 83000000.0; //total no of individuals

int main(){
	double E0 = 30000.0; //initial number of exposed individuals
	double J0 = 9000.0; //initial no of infectious individuals
	double S0 = N - E0 - J0; //initial number of susceptible individuals
	double E;
	double S;
	double J;
	double beta; // equal to R0 * gamma
	double E_new;
	double I_new;
	double S_new;
	double R0[3] = {1.25,1.5,2.0}; //reproduction numbers
	double delta = 0.01;
	double days = 365.0; //number of days 
	double t = 0; //time in days
	FILE* myFile[3] = {fopen("epidemic1.csv", "w") , fopen("epidemic2.csv", "w") , fopen("epidemic3.csv", "w")};

	
	for(int i = 0; i < 3; i++){
		beta = R0[i] * gamma;
		E = E0;
		J = J0;
		S = S0;
		while(t < days){
			S_new = (-beta/N * J * delta + 1) * S;
			E_new = (beta/N * J * S - a * E) * delta + E;
			I_new = (a * E - gamma * J) * delta + J;
			fprintf(myFile[i], "%.10lf \t %.10lf \t %.10lf \t %.10lf  \n", t, S, E, J);
			S = S_new;
			E = E_new;
			J = I_new;
			t += delta;
		}
	fclose(myFile[i]);
	t = 0;
	}
	


	return 0;
}
