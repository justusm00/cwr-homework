#ifndef MYNUMERICS_H
#define MYNUMERICS_H


// Definition of a tuple 
typedef struct {
double x1;
double x2;
}tuple;

/*pointer on function*/
typedef
int ode_func(double, const double[], double[], void * );


// Numerical Integration of scalar 1D function 
double mn_integrate(double left, double right, int N, double integrand(double));



//Implementations of the error function 
double mn_erf_simpson(double x, double delta_x);
double mn_erf_midpoint(double x, double delta_x);



//solutions of a 2nd order polynomial
tuple solve_quadratic(double a, double b , double c);


// derivative
double diff(double x, double delta, double func(double));


//find roots of a polynomial using the Newton-Raphson Algorithm. 
//func: function, x0: initial guess, delta: increment for diff, rel_tol: desired precision, max_iter: stop after this
double find_root(double func(double), double x0, double delta, double rel_tol, int max_iter);

//euler algorithm to solve ode systems
void euler_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params);


//rk2 algorithm to solve ode systems
void rk2_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params);

//rk4 algorithm to solve ode systems
void rk4_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params);


//verlet algorithm to solve ode systems (preserves conserved quantities like energy)
void verlet_step(double t, double delta_t, double y[], ode_func func, int dimension, void* params);


//gaussian elimination
//pivoting is the line of the pivoting element
int gauss(int n, int m, double A[n * m], int pivoting);

//determinant algorithm
double det(int n, double A[n*n], int pivoting);


//system of linear equations
void lgs_solve(int n, double A[n*n], double b[n], double x[n]);

//matrix inverse
void mat_inverse(int n, double A[n*n], double A_inverse[n*n]);

//solve lgs with jacobi method
void lgs_jacobi_solve(int n, double A[], double b[], double x[], int k_max);


//solve lgs gauss seidel method
void lgs_gs_solve(int n, double A[], double b[], double x[], int k_max);

#endif
