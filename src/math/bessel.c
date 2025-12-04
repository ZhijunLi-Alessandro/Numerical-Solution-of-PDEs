/**
 * @file bessel.c
 * @brief Implementations of complex Bessel/Hankel helpers and plane-solution
 *        utilities used by analytic test cases.
 *
 * The low-level functions compute J0 and Y0 via power-series expansions for
 * complex arguments. Higher-level helpers build Hankel combinations and real
 * plane-solution quantities derived from these special functions.
 * 
 * @author Li Zhijun
 * @date 2025-12-02
 */
#include "bessel.h"

/**
 * @brief Compute the Bessel function J_0 for a complex argument using its
 *        power-series expansion.
 *
 * The series used is the standard Bessel J0 Taylor series:
 *   J0(z) = sum_{k=0}^
 *            (-1)^k (z/2)^{2k} / (k!)^2
 * The summation is performed until terms are small relative to the running sum
 * or until a maximum iteration count is reached.
 *
 * @param z Complex argument.
 * @return Complex value of J_0(z).
 */
double complex bessel_J0_complex(double complex z) {
    double complex term = 1.0 + 0.0*I;
    double complex sum = term;

    for (int k = 1; k < 50; k++) {
        term *= -1.0 * (z*z) / (4.0 * k * k);
        sum += term;

        if (cabs(term) < 1e-15 * cabs(sum)) break;
    }
    return sum;
}

/**
 * @brief Compute the Bessel function Y_0 for a complex argument using the
 *        relationship to J_0 and a series with harmonic-number coefficients.
 *
 * The implementation follows the standard expansion for Y0 in terms of J0
 * and a logarithmic/harmonic series. This function assumes `z` is not on the
 * negative real axis branch cut for the complex logarithm. Caller should
 * ensure argument avoids singular points.
 *
 * @param z Complex argument.
 * @return Complex value of Y_0(z).
 */
double complex bessel_Y0_complex(double complex z) {
    const double gamma = 0.5772156649015328606;
    double complex J0 = bessel_J0_complex(z);

    double complex sum = 0.0 + 0.0*I;
    double complex term = 1.0;

    for (int k = 1; k < 50; ++k) {
        term *= -1.0 * (z*z)/(4.0 * k*k);

        double Hk = 0.0;
        for (int m = 1; m <= k; ++m) Hk += 1.0/m;

        sum += term * Hk;

        if (cabs(term) < 1e-15 * cabs(sum)) break;
    }

    return (2.0/M_PI) * ( (gamma + clog(z/2.0)) * J0 - sum );
}

/**
 * @brief Hankel function H_0^{(2)}(z) = J_0(z) - i Y_0(z).
 *
 * @param z Complex argument.
 * @return Complex value of the Hankel function of the second kind.
 */
double complex hankel_H0_2(double complex z) {
    return bessel_J0_complex(z) - I * bessel_Y0_complex(z);
}

/**
 * @brief Evaluate the real plane-solution derived from Hankel function combination.
 *
 * Implementation notes:
 * - The function protects against r==0 by clamping a tiny positive radius.
 * - Uses sqrt(-i) = exp(-i*pi/4) to form the complex argument z = sqrt(-i)*r.
 * - Evaluates H_0^{(2)}(z) and multiplies by exp(i t), returning the real part.
 *
 * @param r Radius (distance from origin). Small values are clamped for stability.
 * @param t Phase parameter (real-valued).
 * @return Real-valued plane-solution (Re{ e^{i t} H_0^{(2)}(sqrt(-i) r) }).
 */
double plane_solution_function(double r, double t) {
    if (r < 1e-8) r = 1e-8;

    /* sqrt(-i) = exp(- i π/4) */
    double complex sqrt_minus_i = cexp(-I * M_PI / 4.0);

    /* argument = sqrt(-i)*r */
    double complex z = sqrt_minus_i * r;

    /* compute hankel_H0_2(z) */
    double complex H0_2 = hankel_H0_2(z);

    /* compute e^{it} */
    double complex eit = cexp(I * t);

    /* value = Re{ e^{it} * H0_2 } */
    double complex val = eit * H0_2;

    return creal(val);
}

/**
 * @brief Approximate average of the plane-solution over a rectangular cell.
 *
 * The formula is derived from asymptotic expansions and includes the
 * Euler--Mascheroni constant. It returns an approximate cell-average used in
 * test RHS/BC assembly.
 *
 * @param hx Cell width in x-direction.
 * @param hy Cell width in y-direction.
 * @param t Phase parameter.
 * @return Approximate average value on the cell.
 */
double average_cell(double hx, double hy, double t) {
    double gamma = 0.5772156649015328606;
    double a = sqrt(hx * hy / M_PI);
    return (cos(t) / 2 + 2 * (log(a / 2) - 0.5 + gamma) * sin(t) / M_PI);
}