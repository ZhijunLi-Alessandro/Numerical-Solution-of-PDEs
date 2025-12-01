#ifndef BESSEL_H
#define BESSEL_H
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// double complex J0_complex(double complex z);
double plane_solution_function(double r, double t);
double average_cell(double hx, double hy, double t);

#endif