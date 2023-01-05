#include <stdio.h>
#include <stdlib.h>
#include <my_numerics.h>






//--------------------MAIN------------------------

int main(void){
	int max_iter = 101;
	double x;
	double x0 = 0.1; //initial value
	double mu[3] = {0.4,0.74,0.77};
	double muuh;
	int iter = 0;
	FILE* myFile[3] = {fopen("pop1.dat", "w") , fopen("pop2.dat", "w") , fopen("pop3.dat", "w")};
	x = x0;
	for(int i = 0; i < 3; i++)
	{
		muuh = mu[i];
		while(iter < max_iter)
		{
			x = 4 * muuh * x * (1 - x);
			fprintf(myFile[i], "%d \t %lf \n", iter,x);
			iter++;
		}
	fclose(myFile[i]);
	iter = 0;
	}
	
	
	



return 0;
}


