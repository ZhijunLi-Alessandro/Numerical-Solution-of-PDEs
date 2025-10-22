/**
 * @file utils.h
 * @brief Header file for utility functions to print vectors and CSR matrices.
 * 
 * This header file declares functions for printing dense vectors,
 * integer vectors, and sparse matrices in Compressed Sparse Row (CSR) format.
 * The printing functions allow for specifying the number of decimal places
 * for floating-point values.
 * @see utils.c
 * @author Li Zhijun
 * @date 2025-10-10
 */
# ifndef UTILS_H
# define UTILS_H
# include "csr.h"

/**
 * @brief Print a dense vector with specified decimal places.
 * @param vec Pointer to the vector.
 * @param n Number of elements in the vector.
 * @param ndec Number of decimal places to print.
 */
void print_vector(const double *vec, int n, int ndec);

/**
 * @brief Print an integer vector.
 * @param vec Pointer to the integer vector.
 * @param n Number of elements in the vector.
 */
void print_int_vector(const int *vec, int n);

void print_matrix(const double **matrix, int rows, int cols, int ndec);

void print_int_matrix(const int **matrix, int rows, int cols);

/**
 * @brief Print a SparseCSR matrix in dense format with specified decimal places.
 * @param matrix Pointer to the SparseCSR matrix.
 * @param ndec Number of decimal places to print.
 * 
 * @note This function is only suitable for CSR matrices with columns sorted within each row.
 */
void print_SparseCSR(const SparseCSR *matrix, int ndec);

/**
 * @brief Print the internal representation of a SparseCSR matrix.
 * @param matrix Pointer to the SparseCSR matrix.
 * @param ndec Number of decimal places to print for values.
 */
void print_SparseCSR_simple(const SparseCSR *matrix, int ndec);

# endif