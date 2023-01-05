#include<stdio.h>
#include<stdlib.h>

// Programm zur numerischen Integration mit der Riemannsumme

double p1(double x);

double integrate(double left, double right, int N, double integrand(double));



int main(){
	FILE* myFile = fopen("data1.txt", "w");
	int N = 10;
	for(int i = 1; i <= 7; i++){
		fprintf(myFile, "%.15lf \n", integrate(0,2,N,p1));
		N *= 10;
	}
	fclose(myFile);
	return 0;
}




double p1(double x){
	double result = x * x * x - x / 2;
	return result;
}

double integrate(double left, double right, int N, double integrand(double)){
	double width = right - left;
	double dx = width / N;
	double result = 0.0;
	for(int i = 0; i < N; i++){
		result += p1(i * dx) * dx;
	}
	return result;
	

}
