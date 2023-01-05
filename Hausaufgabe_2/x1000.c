#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <my_numerics.h>






//--------------------MAIN------------------------

int main(void){
	int max_iter = 1000;
	double x;
	double x0 = 0.01; //initial value
	double mu = 0.4;
	int iter = 0;
	FILE* myFile = fopen("x1000.csv", "w");
	x = x0;
	for(int j = 1; j < max_iter; j++)
	{
		fprintf(myFile,"%lf \t", mu);
		while(iter < max_iter)
		{
			x = 4 * mu * x * (1 - x);
			if(iter > 980)
			{
				fprintf(myFile,"%lf \t", x);
			}
			iter++;
		}
		fprintf(myFile, "\n");
		mu = 0.4 + (1-0.4) * log(j)/log(max_iter);
		iter = 0;
	}
	fclose(myFile);
	
	



return 0;
}


