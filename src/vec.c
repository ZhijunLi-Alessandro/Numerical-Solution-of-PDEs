/**
 * @file vec.c
 * @brief Implementation of basic vector operations.
 * 
 * This file contains the implementation of functions for basic vector operations such as
 * copying, addition, subtraction, scaling, and dot product.
 * @see vec.h
 * @author Li Zhijun
 * @date 2025-10-10
 */
#include "vec.h"

void vec_copy(double *dest, const double *src, int n) {
    for (int i = 0; i < n; i++) {
        dest[i] = src[i];
    }
}

void vec_add(double *a, const double *b, int n) {
    for (int i = 0; i < n; i++) {
        a[i] += b[i];
    }
}

void vec_sub(double *a, const double *b, int n) {
    for (int i = 0; i < n; i++) {
        a[i] -= b[i];
    }
}

void vec_scale(double *a, double scalar, int n) {
    for (int i = 0; i < n; i++) {
        a[i] *= scalar;
    }
}
double vec_dot(const double *a, const double *b, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += a[i] * b[i];
    }
    return result;
}