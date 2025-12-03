#include "bessel.h"

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

double complex hankel_H0_2(double complex z) {
    return bessel_J0_complex(z) - I * bessel_Y0_complex(z);
}

double plane_solution_function(double r, double t) {
    if (r < 1e-8) r = 1e-8;

    // sqrt(-i) = exp(- i π/4)
    double complex sqrt_minus_i = cexp(-I * M_PI / 4.0);

    // argument = sqrt(-i)*r
    double complex z = sqrt_minus_i * r;

    // compute hankel_H0_2(z)
    double complex H0_2 = hankel_H0_2(z);

    // compute e^{it}
    double complex eit = cexp(I * t);

    // f = Re{ e^{it} * J0(z) }
    double complex val = eit * H0_2;

    return creal(val);
}

double average_cell(double hx, double hy, double t) {
    double gamma = 0.5772156649015328606;
    double a = sqrt(hx * hy / M_PI);
    return (cos(t) / 2 + 2 * (log(a / 2) - 0.5 + gamma) * sin(t) / M_PI);
}