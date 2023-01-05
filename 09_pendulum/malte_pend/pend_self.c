/*----------------------------------*\
|  CWR 2021                          |
|  Blatt 5                           |
\*----------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "my_numerics.h"
#include "time.h"
#include "math.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

/*------------------------  PHYSIKALISCHE KONSTANTEN  -------------------------*/

const int N = 1;              // Number of pendulums
const double k = 100;         // Spring constant [N/m]
const double base_length = 1; // zero length of springs [m]
const double mass = 1;        // pendulum mass [kg]

const double v0 = 1;

/*-------------------------  SIMULATIONS-PARAMETER  ---------------------------*/

static const char last_pos_file_name[] = "lastpos.dat";

/*-------------------------  PHYSIKALISCHES SYSTEM  ---------------------------*/

int pendulums_ode(double t, const double y[], double f[], void *params)
{
    // Reset everything
    for (int i = 0; i < 2 * N; ++i)
        f[i] = 0;

    // Positional derivatives: dx/dt = v
    for (int i = 0; i < N; ++i)
        f[i] = y[i + N];

    // Springs between pendulums:
    // Iterate over all pairs
    for (int i = 0; i < N - 1; ++i)
    {
        double delta_x = y[i + 1] - y[i];
        double spring_force = k * (delta_x - base_length);
        f[N + i] += spring_force / mass;
        f[N + i + 1] -= spring_force / mass;
    }

    // Bind first pendulum to origin
    f[N] -= k * y[0] / mass;

    // Return life sign
    return 1;
}

double pendulum_energy(const double y[])
{
    double energy = 0;

    // Kinetic energy
    for (int i = 0; i < N; ++i)
        energy += 0.5 * mass * y[N + i] * y[N + i];

    // Spring energy between pendulums
    for (int i = 0; i < N - 1; ++i)
    {
        double dx = y[i + 1] - y[i] - base_length;
        energy += 0.5 * k * dx * dx;
    }

    // Origin spring
    energy += 0.5 * k * y[0] * y[0];

    return energy;
}

/*-----------------------------------------------------------------------------*/

int main()
{
    /* System vector dimensionality: y = (x1, x2 ... xn, v1, v2 ... vn) */
    int dimension = 2 * N;

    /* ------- Anfangswerte ---------------------------------------------------*/

    /* Allocate system state array and set initial conditions */
    double y_eu[dimension];
    double y_rk2[dimension];

    /* ------- Simulations-Schleife -------------------------------------------*/

    /* Ouput file */
    FILE *last_pos_file_euler = fopen("rk2.dat", "w");
    FILE *last_pos_file_rk2 = fopen("rk4.dat", "w");

    srand(time(NULL));

    /* Calculation quarter period T/4 */
    double T = 2 * M_PI * sqrt(mass / k);
    /* Analytical position after T/4 */
    double x_max = v0 * sqrt(mass / k);

    /* Construct delta_t range */
    int n_exp = 100;                                      // Number of exponents for delta_t
    double min_exp = -8;                                  // Minimal exponent
    double max_exp = 0;                                   // Maximum exponent
    double delta_exp = (max_exp - min_exp) / (n_exp - 1); // Delta between exponents

    /* Goint through all exponents, and to a complete simulation
       for each of the dt's */
    for (int i = 0; i < n_exp; ++i)
    {
        /* Calculate delta_t */
        double current_exp = min_exp + i * delta_exp;
        double delta_t = pow(10, current_exp);

        printf("Calculating now for dt = 10^%.1lf\n", current_exp);

        /* Set initial conditions */
        for (int p = 0; p < N; ++p)
        {
            y_eu[p] = p * base_length;
            y_eu[N + p] = 0;
            y_rk2[p] = p * base_length;
            y_rk2[N + p] = 0;
        }

        y_eu[N] = v0;
        y_rk2[N] = v0;

        /* Do the simulation loop */

        double t = 0;
        double res_eu = 0;
        double res_rk2 = 0;

        double T_max = 1./4 * T;
        int boool = 0;
        while (t <= T_max)
        {

            if (boool) {
                res_eu = (y_eu[0] - x_max) * (y_eu[0] - x_max);
                res_rk2 = (y_rk2[0] - x_max) * (y_rk2[0] - x_max);
                printf("%f, %f\n", T_max, t);
                break;
            }

            if ((T_max - t) < delta_t)
            {
                delta_t = T_max - t;
                boool = 1;
                printf("fg%f, %f\n", T_max, t);
            }
            rk2_step(t, delta_t, y_eu, pendulums_ode, dimension, NULL);
            rk4_step(t, delta_t, y_rk2, pendulums_ode, dimension, NULL);
            t += delta_t;
        }

        /* Output the result */
        fprintf(last_pos_file_euler, "%g\t%g\n", pow(10, current_exp), sqrt(res_eu));
        fprintf(last_pos_file_rk2, "%g\t%g\n", pow(10, current_exp), sqrt(res_rk2));
    }
}
