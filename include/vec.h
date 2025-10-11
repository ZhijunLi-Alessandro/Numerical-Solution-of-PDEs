/**
 * @file vec.h
 * @brief Header file for basic vector operations.
 * 
 * This header file declares functions for basic vector operations such as
 * copying, addition, subtraction, scaling, and dot product.
 * @see vec.c
 * @author Li Zhijun
 * @date 2025-10-10
 */
# ifndef VEC_H
# define VEC_H

/**
 * @brief Copy vector src to dest.
 * @param dest Destination vector.
 * @param src Source vector.
 * @param n Number of elements to copy.
 */
void vec_copy(double *dest, const double *src, int n);

/**
 * @brief Add vector b to vector a (a = a + b).
 * @param a First vector and result vector.
 * @param b Second vector to add.
 * @param n Number of elements in the vectors.
 */
void vec_add(double *a, const double *b, int n);

/**
 * @brief Subtract vector b from vector a (a = a - b).
 * @param a First vector and result vector.
 * @param b Second vector to subtract.
 * @param n Number of elements in the vectors.
 */
void vec_sub(double *a, const double *b, int n);

/**
 * @brief Scale vector a by a scalar (a = a * scalar).
 * @param a Vector to scale.
 * @param scalar Scaling factor.
 * @param n Number of elements in the vector.
 */
void vec_scale(double *a, double scalar, int n);

/**
 * @brief Compute the dot product of vectors a and b.
 * @param a First vector.
 * @param b Second vector.
 * @param n Number of elements in the vectors.
 * @return The dot product (a . b).
 */
double vec_dot(const double *a, const double *b, int n);

# endif