/**
 * @file bessel.h
 * @brief Declarations for Bessel/Hankel-based helper functions used by
 *        analytical / semi-analytical plane solutions in tests and examples.
 *
 * This header exposes routines used to evaluate complex Bessel/Hankel
 * functions (implemented in `src/math/bessel.c`) and higher-level helpers
 * that compute plane-wave-like solutions and cell averages used by the
 * parabolic/heat examples.
 * 
 * @author Li Zhijun
 * @date 2025-12-02
 */
#ifndef BESSEL_H
#define BESSEL_H
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * NOTE: low-level complex Bessel functions are implemented in `bessel.c`.
 * Their signatures are not part of the public API here, but may be
 * re-enabled if needed (commented out below).
 */
// double complex J0_complex(double complex z);

/**
 * @brief Evaluate an analytic plane-like solution built from Hankel functions.
 *
 * The function computes the real part of e^{i t} * H_0^{(2)}(sqrt(-i) * r),
 * where H_0^{(2)} is the Hankel function of the second kind and r is the
 * radial coordinate. The small-radius behaviour is regularized by forcing a
 * minimum r value.
 *
 * @param r Radial coordinate (r >= 0). Small values are clamped for stability.
 * @param t Time-like phase parameter.
 * @return The real-valued plane solution at (r,t).
 */
double plane_solution_function(double r, double t);

/**
 * @brief Approximate average value of the plane solution over a cell.
 *
 * This helper computes an asymptotic average used when writing cell-centered
 * boundary/source terms for tests. The formula uses the Euler--Mascheroni
 * constant and a leading logarithmic dependence on the effective cell area.
 *
 * @param hx Cell width in x-direction.
 * @param hy Cell width in y-direction.
 * @param t Time-like phase parameter.
 * @return Approximate average of the plane solution on a cell of size hx*hy.
 */
double average_cell(double hx, double hy, double t);

#endif