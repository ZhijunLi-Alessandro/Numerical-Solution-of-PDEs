/**
 * @file utils.c
 * @brief Utility functions for printing vectors and SparseCSR matrices.
 * 
 * This file contains functions to print dense vectors, integer vectors,
 * and sparse matrices in Compressed Sparse Row (CSR) format. The printing
 * functions allow for specifying the number of decimal places for floating-point
 * values.
 * 
 * Functions:
 * - print_vector: Print a dense vector with specified decimal places.
 * - print_int_vector: Print an integer vector.
 * - print_SparseCSR: Print a SparseCSR matrix in dense format with specified decimal places.
 * - print_SparseCSR_simple: Print the internal representation of a SparseCSR matrix.
 * 
 * @author Li Zhijun
 * @date 2025-10-10
 */
#include <stdio.h>
#include "utils.h"

void print_vector(const double *vec, int n, int ndec) {
    printf("[");
    for (int i = 0; i < n; i++) {
        printf("%.*f ", ndec, vec[i]);
    }
    printf("]\n");
}

void print_int_vector(const int *vec, int n) {
    printf("[");
    for (int i = 0; i < n; i++) {
        printf("%d ", vec[i]);
    }
    printf("]\n");
}

void print_SparseCSR(const SparseCSR *matrix, int ndec) {
    /*Only Suitable for Columns Sorted CSR Format*/
    for (int i = 0; i < matrix->rows; i++) {
        if (i == 0) {printf("[");}
        else {printf(" ");}
        int k = matrix->row_ptr[i];
        printf("[");
        for (int j = 0; j < matrix->cols; j++) {
            if (k < matrix->row_ptr[i + 1] && matrix->col_ind[k] == j) {
                printf("%.*f ", ndec, matrix->values[k]);
                k++;
            } else {
                printf("%.*f ", ndec, 0.0);
            }
        }
        printf("]");
        if (i == matrix->rows - 1) {printf("]\n");}
        else {printf(",\n");}
    }
}

void print_SparseCSR_simple(const SparseCSR *matrix, int ndec) {
    printf("row_ptr:\n");
    print_int_vector(matrix->row_ptr, matrix->rows + 1);
    printf("col_ind:\n");
    print_int_vector(matrix->col_ind, matrix->nnz);
    printf("values:\n");
    print_vector(matrix->values, matrix->nnz, ndec);
    printf("\n");
}