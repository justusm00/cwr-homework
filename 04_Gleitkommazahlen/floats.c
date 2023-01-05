#include<stdio.h>
#include<stdlib.h>
#include<math.h>


typedef struct {
double x1;
double x2;
}tuple;

tuple solve_quadratic(double a, double b , double c);  //returnt Nullstellen von Polynom 2. Ordnung, nimmt parameter 


//############### MAIN  ######################


int main()
{
	tuple x; //Lösungen
	tuple x_ana; //analytische Lösung
	tuple error; //relative Fehler
	x_ana.x1 = 1e-8;
	x_ana.x2 = 1e8;	
		
	//Parameter für Polynom	
	double a = 1; 
	double b = (-1)*(1e16 + 1) / (1e8);
	double c = 1;
	x = solve_quadratic(a,b,c);
	
	//berechnen der relativen Fehler
	error.x1 = (x_ana.x1 - x.x1) / x_ana.x1; 
	error.x2 = (x_ana.x2 - x.x2) / x_ana.x2;
	
	//ausgabe der Lösungen
	printf("Die Lösungen sind %.15lf und %.15lf \n", x.x1, x.x2);	
	printf("Die relativen Fehler sind %.15lf und %.15lf \n", error.x1, error.x2);
	return 0;
}

//############### solve_quadratic ######################


tuple solve_quadratic(double a, double b , double c)
{	
	tuple x;
	double d = sqrt(b * b - 4 * a * c);	
	double test1;
	double test2;
	
	//Vermeiden, dass Zähler bzw Nenner sehr klein werden
	if(fabsf(-b + d) > 1e-7)
	{
		x.x1 = (-b + d) / (2 * a);
		x.x2 = (2*c) / (-b + d);
	}
	else
	{
		x.x1 = (2*c) / (-b - d);
		x.x2 = (-b - d) / (2 * a);
	}
	test1 = x.x1;
	test2 = x.x2;
	
	//Ausgabe in aufsteigender Reihenfolge
	if(test1 > test2) 
	{
		x.x1 = test2; 
		x.x2 = test1;
	}
	return x;

}
